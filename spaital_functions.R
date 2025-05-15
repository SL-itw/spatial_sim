## covariance matrix ################

c_r_eigen = function(tau, rho, W){


  M = diag(rowSums(W))

  M_vec = diag(M)
  M_in = diag(1/(sqrt(M_vec)))
  matrix_check = M_in %*% W %*% M_in
  eigens = eigen(matrix_check,only.values = T,symmetric = F)$values
  min_eigen = min(eigens)
  max_eigen = max(eigens)

  if (rho >= 1/min_eigen & rho <= 1/max_eigen ){

    Q =  (M - rho*W)
  }else{
    message( paste0("rho should be between ", min_eigen," and ",max_eigen))
    message( paste0("rho pushed to min = ", min_eigen))

    Q =  (M - min_eigen*W)
  }

  # using row-echolon form to find the inverse of a matrix
#matlib::inv(Q) # not all matrixes have an inverse but all have a generalized inverse "Harville - Linear Models and the .. 2028"
matlib::Ginv(tau*Q)
}

############ generating data

gen_data = function(map_data, tau, rho){

  nb <- spdep::poly2nb(map_data, queen = T )


  W = as.matrix(spdep::nb2mat(nb, style = "B",zero.policy = T))

  # removing neighbors with 0 adjacency since they can not contribute and cause singular matrices
 # no_neighborhos = tibble(zero = nb == "0") %>% mutate( id = row_number()) %>% filter(zero == T) %>% pull(id)
  #W_small = W[-no_neighborhos,-no_neighborhos]

  #Q^{-1} = C
  C = c_r_eigen(tau = tau,
                rho = rho,
                W = W
                #W = W_small
                )


  b <- mnormt::rmnorm(n = 1, mean = c(rep(0,dim(W)[1])), varcov = C)
  b
}

## vect version #################

gen_vec = function(map, tau, rho){

  nb <- spdep::poly2nb(map, queen = T )


  W = as.matrix(spdep::nb2mat(nb, style = "B",zero.policy = T))

  # removing neighbors with 0 adjacency since they can not contribute and cause singular matrices
  ids = tibble(zero = nb == "0") %>% mutate( id = row_number())
  neighbors = ids %>% filter(zero == F) %>% pull(id)
  no_neighbors = ids %>% filter(zero == T) %>% pull(id)
  W_small = W[-no_neighbors,-no_neighbors]

  # Q^{-1} = C
  C = c_r_eigen(tau = tau,
                rho = rho,
                #W=W
                W = W_small
                )


  u_i = mnormt::rmnorm(n = 1,
                            mean = c(rep(0,length(neighbors))),
                            varcov = C) %>%
    tibble() %>%
    rename("value" = ".") %>%
    mutate(id = as.numeric(neighbors)) %>%
    right_join(ids %>%
                 mutate(id = as.numeric(id)), by = "id") %>%
    arrange(id)


  u = u_i %>% pull(value)
  v = rnorm(1,0,1)
  b = (sqrt(rho)*u+sqrt(1-rho)*v)/(1/sqrt(tau))

  log_eta = 1 + b

  lambda = exp(log_eta)

  y = rpois(length(lambda),lambda)

  tibble(y = y, lambda = lambda, b = b, u = u, v = v)
}


gen_vec_inla <- function(map, tau, rho) {
  nb <- spdep::poly2nb(map, queen = TRUE)
  W <- as.matrix(spdep::nb2mat(nb, style = "B", zero.policy = TRUE))

  ids <- tibble(zero = apply(W, 1, sum) == 0, id = 1:nrow(W))
  neighbors <- ids %>% filter(!zero) %>% pull(id)
  no_neighbors <- ids %>% filter(zero) %>% pull(id)

  W_small <- if(length(no_neighbors) > 0) W[-no_neighbors, -no_neighbors] else W

  # Sparse Q matrix
  M <- Matrix::Diagonal(x = rowSums(W_small))
  #M <- Matrix::Diagonal(x = rowSums(W))
  Q <- tau*(M - rho * Matrix::Matrix(W_small, sparse = TRUE))
  n <- nrow(W_small)
  I_n <- diag(n)

  # Here, tau_b is the precision (or a hyperparameter for the unstructured part)
  Q11 <- (tau / (1 - rho)) * I_n
  # For Sigma22, note: if you already have Q as tau*(M - rho * W_small), you may adjust it by adding (rho/(1-rho))*I.
  Q22 <- as.matrix(Q) + (rho/(1 - rho)) * I_n
  Q12 <- (-sqrt(tau * rho) / (1 - rho)) * I_n

  # Construct the full block covariance matrix.
  Q_sparse <- rbind(
    cbind(Q11, Q12),
    cbind(Q12, Q22)
  )

  Sigma_sparse = matlib::Ginv(Q_sparse)
  library(MASS)
  # Mean vector of zeros for both parts:
  mu_full <- rep(0, 2 * n)
  joint_sample <- mvrnorm(n = 1, mu = mu_full, Sigma = Sigma_sparse)

  b = joint_sample[1:n]
  u = joint_sample[(n+1):(2*n)]
  # Linear predictor
  eta <- 1 + b
  lambda <- exp(eta)
  # Simulate Poisson counts
  y_obs <- rpois(length(lambda), lambda)
  # Create full vectors for the full map (of length 221)
  total_areas <- nrow(W)
  full_y <- rep(NA, total_areas)
  full_lambda <- rep(NA, total_areas)
  full_b <- rep(NA, total_areas)
  full_u <- rep(NA, total_areas)

  # Fill in the computed values for indices corresponding to "neighbors"
  full_y[neighbors] <- y_obs
  full_lambda[neighbors] <- lambda
  full_b[neighbors] <- b
  full_u[neighbors] <- u

  tibble(y = full_y, lambda = full_lambda, b = full_b, u = full_u)

}

gen_count_data_inla <- function(map, tau, rho) {
  nb <- spdep::poly2nb(map, queen = TRUE)
  W <- as.matrix(spdep::nb2mat(nb, style = "B", zero.policy = TRUE))

  ids <- tibble(zero = apply(W, 1, sum) == 0, id = 1:nrow(W))
  neighbors <- ids %>% filter(!zero) %>% pull(id)
  no_neighbors <- ids %>% filter(zero) %>% pull(id)
  W_small <- if(length(no_neighbors) > 0) W[-no_neighbors, -no_neighbors] else W

  # Sparse Q matrix
   M <- Matrix::Diagonal(x = rowSums(W_small))
  #M <- Matrix::Diagonal(x = rowSums(W))
   Q <- M - rho * Matrix::Matrix(W_small, sparse = TRUE)
  #Q <- M - rho * Matrix::Matrix(W, sparse = TRUE)
  Q_scaled = tau*Q

  # Sample structured spatial effect
  w_1 <- as.vector(inla.qsample(n = 1, Q = Q_scaled ))

  # Create full spatial field with NAs in no-neighbor slots

  # Linear predictor
  eta <- 1+ w_1
  lambda <- exp(eta)

  # Simulate Poisson counts
  y_obs <- rpois(length(lambda), lambda)

  # Create full vectors for the full map (of length 221)
  total_areas <- nrow(W)
  full_y <- rep(NA, total_areas)
  full_lambda <- rep(NA, total_areas)
  full_b <- rep(NA, total_areas)

  # Fill in the computed values for indices corresponding to "neighbors"
  full_y[neighbors] <- y_obs
  full_lambda[neighbors] <- lambda
  full_b[neighbors] <- w_1

  tibble(y = full_y, lambda = full_lambda, b = full_b)
}


gen_count_data_inla <- function(map, tau, rho, beta_0 = 1, beta = 0, covariate = NULL,
                                neighbor_type = "queen", neighbor_param = 1) {

  if(neighbor_type == "queen") {
    # Get first order neighbors using queen contiguity:
    nb <- poly2nb(map, queen = TRUE)
    # If a higher order is desired (e.g., 2 or 3), union the neighbors at each lag:
    if(neighbor_param > 1) {
      nb_list <- nb
      for(j in 2:neighbor_param) {
        nb_j <- nblag(nb, j)[[j]]
        nb_list <- union.nb(nb_list, nb_j)
      }
      nb <- nb_list
    }
  } else if(neighbor_type == "distance") {
    # For a distance-based matrix, assume the map has coordinates.
    # (If your map is an sf object in longlat, set longlat = TRUE)
    coords <- sp::coordinates(map)
    nb <- dnearneigh(coords, 0, neighbor_param, longlat = FALSE)
  } else {
    stop("neighbor_type must be either 'queen' or 'distance'.")
  }

  # Convert neighbors to a binary weight matrix:
  W <- as.matrix(nb2mat(nb, style = "B", zero.policy = TRUE))

  # Identify areas with no neighbors:
  ids <- tibble(zero = apply(W, 1, sum) == 0, id = 1:nrow(W))
  neighbors <- ids %>% filter(!zero) %>% pull(id)
  no_neighbors <- ids %>% filter(zero) %>% pull(id)
  W_small <- if(length(no_neighbors) > 0) W[-no_neighbors, -no_neighbors] else W
# check if generalized inverse can handle 0 neighbors
  # Construct the precision matrix:
  M <- Matrix::Diagonal(x = rowSums(W_small))
  Q <- M - rho * Matrix::Matrix(W_small, sparse = TRUE)


  # Sample the structured spatial effect (φ)
  w <- as.vector(inla.qsample(n = 1, Q = Q_scaled))
  # check if this compares rmvn is comparable to inla.qsample

  # Create a full spatial field: fill in NA for areas with no neighbors.
  w_1 <- rep(NA, nrow(W))
  w_1[neighbors] <- w

  # Create covariate if not provided (intercept only)
  if(is.null(covariate)) {
    covariate <- rep(0, length(w_1))
  }

  # Compute the linear predictor and intensity:
  eta <- beta_0 + beta * covariate + w_1
  lambda <- exp(eta)

  # Simulate Poisson counts for Y:
  y_obs <- rpois(length(lambda), lambda)

  # Return a tibble with an added region ID:
  tibble(region = 1:nrow(W), y = y_obs, lambda = lambda, w_1 = w_1)
}

# Simulation wrapper: run the function n_sim times and collect the results
simulate_inla <- function(map, tau, rho, beta_0 = 1, beta = 0, covariate = NULL,
                          neighbor_type = "queen", neighbor_param = 1, n_sim = 100) {
  results_list <- vector("list", n_sim)
  for(i in 1:n_sim) {
    results_list[[i]] <- gen_count_data_inla(map, tau, rho, beta_0, beta, covariate,
                                             neighbor_type, neighbor_param) %>%
      mutate(sim = i)  # add simulation number
  }
  all_results <- bind_rows(results_list)
  return(all_results)
}

# Summarize the simulation results:
summarize_simulation <- function(sim_results) {
  sim_results %>%
    group_by(region) %>%
    summarise(
      y_median = median(y, na.rm = TRUE),
      y_lower  = quantile(y, 0.025, na.rm = TRUE),
      y_upper  = quantile(y, 0.975, na.rm = TRUE),
      w_1_median = median(w_1, na.rm = TRUE),
      w_1_lower  = quantile(w_1, 0.025, na.rm = TRUE),
      w_1_upper  = quantile(w_1, 0.975, na.rm = TRUE)
    )
}

# Sim functions ----------------

## Adjacentcy matrices

get_adj <- function(map,
                         style = c("queen", "distance"),
                         order = 1,
                         dist_threshold = NULL) {
  # build nb based on style
  if (style == "queen") {
    nb1 <- poly2nb(map, queen = TRUE)
    # for higher‐order contiguity use nblag
    if (order > 1) {
      nb_list <- nblag(nb1, maxlag = order)
      nb_use  <- nblag_cumul(nb_list)
    } else {
      nb_use <- nb1
    }
  } else if (style == "distance") {
    # need centroids and a numeric threshold in map's units
    pts <- st_centroid(map)
    coords <- st_coordinates(pts)
    nb_use <- dnearneigh(coords, d1 = 0, d2 = dist_threshold, longlat = FALSE)
  } else {
    stop("style must be 'queen' or 'distance'")
  }

  # convert to adjacency and drop islands
  W      <- nb2mat(nb_use, style = "B", zero.policy = TRUE) %>% as.matrix()
  zero_i <- which(rowSums(W) == 0)
  if (length(zero_i) > 0) {
    W_small <- W[-zero_i, -zero_i]
  } else {
    W_small <- W
  }

  keep =setdiff(seq_len(nrow(map)), zero_i)
  list(W=W, W_small = W_small, nb_use = nb_use, zero_i = zero_i, keep = keep)
}

scale_q <- function(W_small, mc = 200) {
  n_s   <- nrow(W_small)
  M   <- diag(rowSums(W_small))
  Q   <- M - W_small + diag(1e-6, n_s)    # jitter for stability
  R   <- chol(Q)
  # trace‐estimate:
  trQinv <- mean(replicate(mc, {
    z <- rnorm(n_s)
    fs <- forwardsolve(t(R), z)
    bs <- backsolve(R, fs)
    sum(z*bs)
  }))
  c     = trQinv/ n_s
  Q_s = c * Q
  list(Q_s = Q_s, n_s = n_s)
}

get_cholesky_bym2 <- function(map,
                              style,
                              order ,
                              dist_threshold ,
                              tau ,
                              rho ) {
  # adjacency
  adjacency_list = get_adj(map = map,
                           style = style,
                           order = order,
                           dist_threshold = dist_threshold)

  W_small = adjacency_list[["W_small"]]

  # precision blocks
  p_data = scale_q(W_small = W_small)
  Q_u_sc = p_data[["Q_s"]]
  n_s = p_data[["n_s"]]
  I_n = diag(n_s)

  Q11      <- (tau / (1 - rho)) * I_n
  Q22      <- Q_u_sc + (rho / (1 - rho)) * I_n
  Q12      <- (-sqrt(tau * rho) / (1 - rho)) * I_n

  Q_block  <- rbind(
    cbind(Q11, Q12),
    cbind(Q12, Q22)
  )

  R_t = t(chol(Q_block + diag(1e-6, 2*n_s)))
  # returns transpose for forward substitution
  list(R_t = R_t,
       W = adjacency_list[["W"]],
       W_small = W_small,
       nb_use = adjacency_list[["nb_use"]],
       zero_i = adjacency_list[["zero_i"]],
       n_s = n_s,
       keep = adjacency_list[["keep"]])
}



get_b = function(map,
                 style,
                 order ,
                 dist_threshold ,
                 tau ,
                 rho){

 data = get_cholesky_bym2(map = map,
                   style = style,
                   order = order,
                   dist_threshold = dist_threshold,
                   tau = tau,
                   rho = rho)

  R_t <- data[["R_t"]]
  n_s <- data[["n_s"]]
  keep <- data[["keep"]]
  # 3. simulate
  rn       <- rnorm(2 * n_s)
  c_vec    <- forwardsolve(R_t, rn)     # solves L · c_vec = rn
  b_small  <- c_vec[1:n_s]

  # 4. “re-inflate” back to length nrow(map)
  b_full <- rep(NA_real_, nrow(map))
  b_full[keep] <- b_small

  b_full
}

## get BYM2 model Q





compute_bym2_chol <- function(W_small, rho, tau) {
  n   <- nrow(W_small)
  I_n <- diag(n)
  # scaled unstructured precision
  Q_u_sc <- tau * (diag(rowSums(W_small)) - rho * W_small)
  Q11 <- (tau      / (1 - rho)) * I_n
  Q22 <- Q_u_sc + (rho / (1 - rho)) * I_n
  Q12 <- (-sqrt(tau * rho) / (1 - rho)) * I_n
  Q_block <- rbind(
    cbind(Q11, Q12),
    cbind(Q12, Q22)
  )
  t(chol(Q_block))
}

### parsing function

# helper to parse the name into its components
parse_spec <- function(spec_name) {
  # 1) the map‐prefix (e.g. “bronx” from “bronx_d1mi_r0.9_t1.5_map_bx”)
  map_name <- str_remove(spec_name, "_(?:q\\d|d[0-9.]+mi)_.*$")

  # 2) the neighbour chunk (q1, q2, q3 or d0.5mi, d1.0mi, etc.)
  nb_chunk <- str_extract(spec_name, "(?:q\\d|d[0-9.]+mi)")

  # 3) style / order / dist_threshold
  if (str_starts(nb_chunk, "q")) {
    style          <- "queen"
    order          <- as.integer(str_remove(nb_chunk, "^q"))
    dist_threshold <- NA_real_
  } else {
    style          <- "distance"
    order          <- NA_integer_
    # extract the number between "d" and "mi"
    dist_threshold <- as.numeric(str_extract(nb_chunk, "(?<=d)[0-9.]+(?=mi)"))
  }

  # 4) rho & tau via look‐around so we never grab the "_map_" part
  rho <- as.numeric( str_extract(spec_name, "(?<=_r)[0-9.]+(?=_)") )
  tau <- as.numeric( str_extract(spec_name, "(?<=_t)[0-9.]+(?=_map_|$)") )

  # 5) the map file suffix (map_m, map_bx, etc.)
  map_file_name <- str_extract(spec_name, "map_[^_]+$")

  list(
    map_name       = map_name,
    style          = style,
    order          = order,
    dist_threshold = dist_threshold,
    rho            = rho,
    tau            = tau,
    map_file_name  = map_file_name
  )
}


# simulation + summarization for a single R_t
simulate_one <- function(spec_name) {

  ps = parse_spec(spec_name)
  map_file_name = ps$map_file_name
  map = qs::qread(paste0("./qs_data/map_sf/",map_file_name,".qs"))
  pop_full <- map %>% st_drop_geometry() %>% pull(pop)
  n_map    <- length(pop_full)

  if(ps$style == "queen"){
    W = qs::qread(paste0("./qs_data/qs_adj/",map_file_name,"_",ps$style,"_o",ps$order,".qs"))$W_small}else {
      W = qs::qread(paste0("./qs_data/qs_adj/",map_file_name,"_",ps$style,"_d",ps$dist_threshold,"mi", ".qs"))$W_small
    }

  keep <- which(rowSums(W) != 0)
  n_s  <- length(keep)
  pop_s <- pop_full[keep]

  # write an INLA graph file and read it back in
  spec_g = paste0("./graphs/ ",ps$map_name,"_",ps$style,"_",ps$order,"_",ps$dist_threshold,".graph")
  graph <- inla.read.graph(spec_g)
  # formula for data check
  formula <- y ~ 1 +
    f(idx,
      model       = "bym2",
      graph       = graph,
      scale.model = TRUE,
      hyper = list(
        prec = list(
          # Log‐Gamma(1,0.001) is very flat ⇒ vague prior on τ
          prior = "loggamma",
          param = c(1, 0.001)
        ),
        phi  = list(prior = "logitbeta", param = c(1, 1))
      )
    )

  # prepare storage
  b_store      <- matrix(NA, nrow = n_map, ncol = r_iter)
  lambda_store <- matrix(NA, nrow = n_map, ncol = r_iter)
  y_store      <- matrix(NA, nrow = n_map, ncol = r_iter)

  rho_store      <- matrix(NA, nrow = n_map, ncol = r_iter)
  tau_store <- matrix(NA, nrow = n_map, ncol = r_iter)
  lambda_est_store      <- matrix(NA, nrow = n_map, ncol = r_iter)
  yhat_store      <- matrix(NA, nrow = n_map, ncol = r_iter)


    chol_name = paste0("./qs_data/qs_chol/",ps$map_name,"_",ps$style,"_",ps$order,"_",ps$dist_threshold,"_",ps$rho,"_",ps$tau,".qs")


  R_t = qs::qread(paste0(chol_name))

  rho_iter <- numeric(r_iter)
  tau_iter <- numeric(r_iter)


  # simulate
  for (i in seq_len(r_iter)) {
    # draw from the Gaussian with precision Q = L Lᵀ, where L = R_t
    rn <- rnorm(2*n_s)
    c_vec <- forwardsolve(R_t, rn)         # solves L · c_vec = rn
    b   <- c_vec[1:n_s]
    #log(pop_s/1000)
    # linear predictor + Poisson
    z  <-   B_0 + b
    lambda <- exp(z)
    y      <- rpois(n_s, lambda)



    # expand to full map
    b_full      <- rep(NA, n_map);      b_full[keep]      <- b
    lambda_full <- rep(NA, n_map); lambda_full[keep] <- lambda
    y_full      <- rep(NA, n_map);      y_full[keep]      <- y

    #estimate prameters
    dat <- tibble(idx = seq_len(n_s), y = y)
    res <- inla(
      formula,
      family            = "poisson",
      data              = dat,
      control.predictor = list(compute = TRUE),
      control.compute   = list(dic = TRUE, cpo = TRUE)
    )


    # extract scalar hyperpars
    hyp    <- res$summary.hyperpar
    i_rho <- grep("phi",  rownames(hyp), ignore.case = TRUE)[1]
    i_tau <- grep("prec", rownames(hyp), ignore.case = TRUE)[1]



      rho_iter[i] <- hyp[i_rho, "0.5quant"]
      tau_iter[i] <- hyp[i_tau, "0.5quant"]


    lambda_est_full <-  rep(NA, n_map); lambda_est_full[keep] <- exp(res$summary.linear.predictor$`0.5quant`)
    yhat_full      <- rep(NA, n_map);      yhat_full[keep]      <- res$summary.fitted.values$`0.5quant`

    b_store[,i]      <- b_full
    lambda_store[,i] <- lambda_full
    y_store[,i]      <- y_full

    lambda_est_store[,i]      <- lambda_est_full
    yhat_store[,i]      <- yhat_full
  }


## summarization helper
summarize_mat <- function(mat, prefix) {
    tibble(
      median = apply(mat, 1, median,    na.rm = TRUE),
      q025   = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
      q975   = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
    ) %>%
      rename_with(~ paste0(prefix, "_", .), everything())
  }

  b_sum      <- summarize_mat(b_store,      "b")
  lambda_est_sum <- summarize_mat(lambda_est_store, "lambda_hat")
  lambda_sum <- summarize_mat(lambda_store, "lambda")
  y_sum      <- summarize_mat(y_store,      "y")
  yhat_sum      <- summarize_mat(yhat_store,      "yhat")

  # bind together, add identifiers
  tibble(
    map_spec = spec_name,
    area_id  = seq_len(n_map),
    rho_median =  median(rho_iter),
    rho_q025   = quantile(rho_iter, .025),
    rho_q975   = quantile(rho_iter, .975),
    rho_bias = (sum(  rho_iter - as.numeric(ps$rho) ))/(r_iter),
    rho_mse = (sum( (( rho_iter - as.numeric(ps$rho) )^2 )))/(r_iter-1),
    tau_median = median(tau_iter),
    tau_q025   = quantile(tau_iter, .025),
    tau_q975   = quantile(tau_iter, .975),
    tau_bias = (sum(  tau_iter - as.numeric(ps$tau) ))/(r_iter),
    tau_mse = (sum( (( tau_iter - as.numeric(ps$tau) )^2 )))/(r_iter-1)
  ) %>%
    bind_cols(lambda_sum,
              lambda_est_sum,y_sum, yhat_sum, b_sum)
}


# simulation + summarization for a single R_t
simulate_pc <- function(spec_name) {

  ps = parse_spec(spec_name)
  map_file_name = ps$map_file_name
  map = qs::qread(paste0("./qs_data/map_sf/",map_file_name,".qs"))
  pop_full <- map %>% st_drop_geometry() %>% pull(pop)
  n_map    <- length(pop_full)

  if(ps$style == "queen"){
    W = qs::qread(paste0("./qs_data/qs_adj/",map_file_name,"_",ps$style,"_o",ps$order,".qs"))$W_small}else {
      W = qs::qread(paste0("./qs_data/qs_adj/",map_file_name,"_",ps$style,"_d",ps$dist_threshold,"mi", ".qs"))$W_small
    }

  keep <- which(rowSums(W) != 0)
  n_s  <- length(keep)
  pop_s <- pop_full[keep]

  # write an INLA graph file and read it back in
  spec_g = paste0("./graphs/ ",ps$map_name,"_",ps$style,"_",ps$order,"_",ps$dist_threshold,".graph")
  graph <- inla.read.graph(spec_g)
  # formula for data check
  formula <- y ~ 1 +
    f(idx,
      model       = "bym2",
      graph       = graph,
      scale.model = TRUE,
      hyper = list(
        prec = list(prior = "pc.prec", param = c(1, 0.9)),     # P(σ > 1) = 0.01
        phi  = list(prior = "pc",      param = c(0.5, 0.1))     # P(φ < 0.5) = 2/3
      )
    )

  # prepare storage
  b_store      <- matrix(NA, nrow = n_map, ncol = r_iter)
  lambda_store <- matrix(NA, nrow = n_map, ncol = r_iter)
  y_store      <- matrix(NA, nrow = n_map, ncol = r_iter)

  rho_store      <- matrix(NA, nrow = n_map, ncol = r_iter)
  tau_store <- matrix(NA, nrow = n_map, ncol = r_iter)
  lambda_est_store      <- matrix(NA, nrow = n_map, ncol = r_iter)
  yhat_store      <- matrix(NA, nrow = n_map, ncol = r_iter)


  chol_name = paste0("./qs_data/qs_chol/",ps$map_name,"_",ps$style,"_",ps$order,"_",ps$dist_threshold,"_",ps$rho,"_",ps$tau,".qs")


  R_t = qs::qread(paste0(chol_name))

  rho_iter <- numeric(r_iter)
  tau_iter <- numeric(r_iter)


  # simulate
  for (i in seq_len(r_iter)) {
    # draw from the Gaussian with precision Q = L Lᵀ, where L = R_t
    rn <- rnorm(2*n_s)
    c_vec <- forwardsolve(R_t, rn)         # solves L · c_vec = rn
    b   <- c_vec[1:n_s]
    u   <- c_vec[(n_s+1):(n_s*2)]
    #log(pop_s/1000)
    # linear predictor + Poisson
    z  <-   B_0 + b
    lambda <- exp(z)
    y      <- rpois(n_s, lambda)



    # expand to full map
    b_full      <- rep(NA, n_map);      b_full[keep]      <- b
    lambda_full <- rep(NA, n_map); lambda_full[keep] <- lambda
    y_full      <- rep(NA, n_map);      y_full[keep]      <- y

    #estimate prameters
    dat <- tibble(idx = seq_len(n_s), y = y)
    res <- inla(
      formula,
      family            = "poisson",
      data              = dat,
      control.predictor = list(compute = TRUE),
      control.compute   = list(dic = TRUE, cpo = TRUE, config = T)
    )

 pred = inla.posterior.sample(result = res,n = 200 )
 lam.rep <- vector("list", 200)
 rho.rep <- numeric(200)
 tau.rep <- numeric(200)

 for(m in seq_len(200)) {
   h    <- pred[[m]]$hyperpar
   ##    - INLA stores log(τ) under something with "Prec" in the name,
   ##      and logit(φ) under something with "Phi" in the name
   lp   <-   h[ grep("prec", names(h), ignore.case=TRUE)[1] ]
   lphi <-   h[ grep("phi",  names(h), ignore.case=TRUE)[1] ]

   ## back‐transform
   tau.rep[m] <- exp(lp)       # τ = exp(log(τ))
   rho.rep[m] <- plogis(lphi)  # φ = logistic(logit(φ))

   samp   <- pred[[m]]$latent
   ηm     <- samp[grep("^Predictor", rownames(samp))]
   λm     <- exp(ηm)
   lam.rep[[m]] <- λm# posterior‐draw of the Poisson means
   y.rep[[m]] <- rpois(length(λm), λm)
 }

 wide.pred <- bind_cols(
   dat,
   do.call(cbind,lam.rep)
 )


    # extract scalar hresult = # extract scalar hyperpars
    hyp    <- res$summary.hyperpar
    i_rho <- grep("phi",  rownames(hyp), ignore.case = TRUE)[1]
    i_tau <- grep("prec", rownames(hyp), ignore.case = TRUE)[1]



    rho_iter[i] <- hyp[i_rho, "0.5quant"]
    tau_iter[i] <- hyp[i_tau, "0.5quant"]


    lambda_est_full <-  rep(NA, n_map); lambda_est_full[keep] <- exp(res$summary.linear.predictor$`0.5quant`)
    yhat_full      <- rep(NA, n_map);      yhat_full[keep]      <- res$summary.fitted.values$`0.5quant`

    b_store[,i]      <- b_full
    lambda_store[,i] <- lambda_full
    y_store[,i]      <- y_full

    lambda_est_store[,i]      <- lambda_est_full
    yhat_store[,i]      <- yhat_full
  }


  ## summarization helper
  summarize_mat <- function(mat, prefix) {
    tibble(
      median = apply(mat, 1, median,    na.rm = TRUE),
      q025   = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
      q975   = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
    ) %>%
      rename_with(~ paste0(prefix, "_", .), everything())
  }

  b_sum      <- summarize_mat(b_store,      "b")
  lambda_est_sum <- summarize_mat(lambda_est_store, "lambda_hat")
  lambda_sum <- summarize_mat(lambda_store, "lambda")
  y_sum      <- summarize_mat(y_store,      "y")
  yhat_sum      <- summarize_mat(yhat_store,      "yhat")

  # bind together, add identifiers
  tibble(
    map_spec = spec_name,
    area_id  = seq_len(n_map),
    rho_median =  median(rho_iter),
    rho_q025   = quantile(rho_iter, .025),
    rho_q975   = quantile(rho_iter, .975),
    rho_bias = (sum(  rho_iter - as.numeric(ps$rho) ))/(r_iter),
    rho_mse = (sum( (( rho_iter - as.numeric(ps$rho) )^2 )))/(r_iter-1),
    tau_median = median(tau_iter),
    tau_q025   = quantile(tau_iter, .025),
    tau_q975   = quantile(tau_iter, .975),
    tau_bias = (sum(  tau_iter - as.numeric(ps$tau) ))/(r_iter),
    tau_mse = (sum( (( tau_iter - as.numeric(ps$tau) )^2 )))/(r_iter-1)
  ) %>%
    bind_cols(lambda_sum,
              lambda_est_sum,y_sum, yhat_sum, b_sum)
}


simulate_moran <- function(rho, tau = 1.5, lw) {
  # build block precision
  Q11 <- (tau / (1 - rho)) * I_n
  Q22 <- Q_star+ (rho / (1 - rho)) * I_n
  Q12 <- (-sqrt(tau * rho) / (1 - rho)) * I_n

  Q_block <- rbind(
    cbind(Q11, Q12),
    cbind(Q12, Q22)
  )+Diagonal(2*n_ws, 1e-8)

  # add tiny ridge for numerical stability
  R <- Matrix::chol(Q_block, pivot = FALSE)
  # simulate joint random vector
  z <- rnorm(2 * n_ws)
  # solve Lᵀ · w = z
  w <- backsolve(R,z, transpose = T)

  # extract latent b
  b    <- w[1:n_ws]
  y    <- rpois(n_ws, lambda = exp(1 + b))

  # Moran's I tests
  mb   <- moran.test(b,  lw)$estimate[["Moran I statistic"]]
  my   <- moran.test(y,  lw)$estimate[["Moran I statistic"]]

  c(moran_b = mb, moran_y = my)
}

simulate_moran_b <- function(rho, tau = 1.5, lw) {

  Prec_b  <- tau * ((1 -rho) * diag(n_ws) +rho * Q_star)
  Sigma_b <- solve(Prec_b)
  b <- as.numeric(mvtnorm::rmvnorm(1, sigma = Sigma_b))
  y    <- rpois(n_ws, lambda = exp(1 + b))

  # Moran's I tests
  mb   <- moran.test(b,  lw)$estimate[["Moran I statistic"]]
  my   <- moran.test(y,  lw)$estimate[["Moran I statistic"]]

  list(moran_b = mb, moran_y = my)
}

simulate_moran_a <- function(rho, tau = 1.5, lw) {

  Prec_b  <- tau * ((1 -rho) * diag(n_ws) +rho * Q_star)
  b <- forwardsolve(t(chol(Prec_b)),rnorm(nrow(Prec_b)))
  y    <- rpois(n_ws, lambda = exp(1 + b))

  # Moran's I tests
  mb   <- moran.test(b,  lw)$estimate[["Moran I statistic"]]
  my   <- moran.test(y,  lw)$estimate[["Moran I statistic"]]

  list(moran_b = mb, moran_y = my)
}

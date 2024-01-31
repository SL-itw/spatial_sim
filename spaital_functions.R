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

    E_inv =  M - rho*W
  }else{
    message( paste0("rho should be between ", min_eigen," and ",max_eigen))
    message( paste0("rho pushed to min = ", min_eigen))

    E_inv =  M - min_eigen*W
  }



  tau^2 * E_inv # not all matrixes have an inverse but all have a generalized inverse "Harville - Linear Models and the .. 2028"

}

############ generating data

gen_data = function(map_data, tau, rho){

  nb <- spdep::poly2nb(map_data, queen = T )


  W = as.matrix(spdep::nb2mat(nb, style = "B",zero.policy = T))

  # removing neighbors with 0 adjacency since they can not contribute and cause singular matrices
  no_neighborhos = tibble(zero = nb == "0") %>% mutate( id = row_number()) %>% filter(zero == T) %>% pull(id)
  W_small = W[-no_neighborhos,-no_neighborhos]

  C = c_r_eigen(tau = tau,
                rho = rho,
                W = W_small)


  gen_data <- mnormt::rmnorm(n = 10, mean = c(rep(0,dim(W)[1])), varcov = C)
  gen_data
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

  C = c_r_eigen(tau = tau,
                rho = rho,
                W = W_small)


  gen_data = mnormt::rmnorm(n = 1,
                            mean = c(rep(0,length(neighbors))),
                            varcov = C) %>%
    tibble() %>%
    rename("value" = ".") %>%
    mutate(id = as.numeric(neighbors)) %>%
    right_join(ids %>%
                 mutate(id = as.numeric(id)), by = "id") %>%
    arrange(id)


  obs_vec = gen_data %>% pull(value)

  obs_vec
}


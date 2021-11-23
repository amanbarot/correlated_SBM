library(R.utils)
library(NlcOptim)
library(assertthat)
library(caret)
library(pROC)
library(MCMCpack)

pow_law_sample <- function(n, x_min, x_max, eta){
  #n, number of samples
  #eta is the exponent of the power law distribution, say eta = 2
  #xmin and xmax control the min and max degree. If we are normalizing the samples, this controls the range of the norm. samples.
  u <- runif(n)
  my_samples <- (x_min^(1-eta) - u*(x_min^(1-eta) - x_max^(1-eta)))^(1/(1-eta))
  return(n * my_samples/sum(my_samples))
}

corr_dcsbm_latent <- function(p1, p2, rho, deg1, deg2, blocks){
  #assumes parameters provided are feasible
  n <- length(blocks)
  K <- max(blocks)
  deg_prod1 <- deg1 %*% t(deg1)
  diag(deg_prod1) <- 0
  deg_prod2 <- deg2 %*% t(deg2)
  diag(deg_prod2) <- 0
  A1 <- matrix(0, nrow = n, ncol = n)
  A2 <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if(i >= j) next
      #browser()
      temp_par <- rep(0, 5)
      #order is 10, 01, 11, 00, and the last term is dev. from independence
      if(rho[blocks[i], blocks[j]] >= 0){
        temp_par[5] <- rho[blocks[i], blocks[j]] * (min(p1[blocks[i], blocks[j]] * deg_prod1[i,j], p2[blocks[i], blocks[j]] * deg_prod2[i,j])) -
          rho[blocks[i], blocks[j]] * p1[blocks[i], blocks[j]] * deg_prod1[i,j] * p2[blocks[i], blocks[j]] * deg_prod2[i,j]
      }
      else{
        temp_par[5] <- -1 * abs(rho[blocks[i], blocks[j]]) * p1[blocks[i], blocks[j]] * deg_prod1[i,j] * p2[blocks[i], blocks[j]] * deg_prod2[i,j]
        if(p1[blocks[i], blocks[j]] * deg_prod1[i,j] + p2[blocks[i], blocks[j]] * deg_prod2[i,j] - 1 >= 0){
          temp_par[5] <- temp_par[5] + abs(rho[blocks[i], blocks[j]]) * 
            (p1[blocks[i], blocks[j]] * deg_prod1[i,j] + p2[blocks[i], blocks[j]] * deg_prod2[i,j] - 1)
        }
      }
      temp_par[1] <- p1[blocks[i], blocks[j]] * deg_prod1[i,j] *
        (1 - p2[blocks[i], blocks[j]] * deg_prod2[i,j]) -
        temp_par[5]
      temp_par[2] <- p2[blocks[i], blocks[j]] * deg_prod2[i,j] *
        (1 - p1[blocks[i], blocks[j]] * deg_prod1[i,j]) -
        temp_par[5]
      temp_par[3] <- p1[blocks[i], blocks[j]] * deg_prod1[i,j] * p2[blocks[i], blocks[j]] * deg_prod2[i,j] +
        temp_par[5]
      temp_par[4] <- (1 - p1[blocks[i], blocks[j]] * deg_prod1[i,j]) *
        (1 - p2[blocks[i], blocks[j]] * deg_prod2[i,j]) +
        temp_par[5]
      #parameters defined. now sample the edge
      temp_par <- pmin(temp_par, 1) #to prevent rounding errors
      temp_par <- pmax(temp_par, 0) 
      if(sum(temp_par[1:4] < 0) > 0){ #a check. perhaps remove later.
        browser()
      }
      temp_edg <- c(rmultinom(1, 1, temp_par[1:4])) 
      #order is 10,01,11,00    
      if(temp_edg[1] == 1){
        A1[i,j] = 1
        A1[j,i] = 1
      }
      
      if (temp_edg[2] == 1) {
        #browser()
        A2[i,j] = 1
        A2[j,i] = 1
      }
      
      if(temp_edg[3] == 1){
        A1[i,j] = 1
        A2[i,j] = 1
        A1[j,i] = 1
        A2[j,i] = 1
      }
    }
  }
  assert_that(sum(abs(A1 - t(A1))) == 0)
  assert_that(sum(abs(A2 - t(A2))) == 0)
  return(list(A1, A2, deg1, deg2))
}

edge_idx <- function(edg, n){
  row <- 1 + ((edg - 1) %% n)
  col <- 1 + ((edg - 1) %/% n)
  return(c(row, col))
}

check_feas_par <- function(p1, p2, rho, deg1, deg2, blocks){
  #checks if p1, p2 and rho are valid parameters given the degrees and blocks
  K <- max(blocks)
  n <- length(blocks)
  deg_prod1 <- deg1 %*% t(deg1)
  diag(deg_prod1) <- 0
  deg_prod2 <- deg2 %*% t(deg2)
  diag(deg_prod2) <- 0
  flag <- T #indicates whether constraints are violated
  for (i in 1:K) {
    for (j in 1:K) {
      if(p1[i,j] < 0 | p2[i,j] < 0 | rho[i,j] < -1){
        print("Lower bound violated")
        flag <- F
      }
      if(p1[i,j] > 1/max(deg_prod1[blocks == i, blocks == j]) | p2[i,j] > 1/max(deg_prod2[blocks == i, blocks == j])){
        print("p1 or p2 violate upper bound")
        flag <- F
      }
      if(rho[i,j] > 1){
        print("Upper bound violated")
        flag <- F
      }
    }  
  }
  return(flag)
}

init_pars_latent <- function(deg1, deg2, blocks, init_rho = NULL){
  n <- length(deg1)
  K <- max(blocks)
  init_p1 <- matrix(0, nrow = K, ncol = K)
  init_p2 <- matrix(0, nrow = K, ncol = K)
  deg_prod1 <- deg1 %*% t(deg1)
  diag(deg_prod1) <- 0
  deg_prod2 <- deg2 %*% t(deg2)
  diag(deg_prod2) <- 0
  if(is.null(init_rho)){
    init_rho <- matrix(0.5, nrow = K, ncol = K)
  }
  for (i in 1:K) {
    for (j in 1:K) {
      init_p1[i,j] <- 0.5 * 1/max(deg_prod1[blocks == i, blocks == j])
      init_p2[i,j] <- 0.5 * 1/max(deg_prod2[blocks == i, blocks == j])
    }
  }
  #browser()
  return(list(init_p1, init_p2, init_rho))
}

estim_dcsbm_latent <- function(A1, A2, deg1, deg2, blocks, fold = NULL, estim_fn = block_pair_estim_latent, lam=NULL){
  n <- dim(A1)[1]
  K <- max(blocks)
  if(is.null(fold)) fold <- which(upper.tri(A1))
  fold_mat <- matrix(F, nrow = n, ncol = n)
  fold_mat[fold] <- T
  deg_prod1 <- deg1 %*% t(deg1)
  diag(deg_prod1) <- 0
  deg_prod2 <- deg2 %*% t(deg2)
  diag(deg_prod2) <- 0
  sol_p1 <- matrix(0, nrow = K, ncol = K)
  sol_p2 <- matrix(0, nrow = K, ncol = K)
  sol_rho <- matrix(0, nrow = K, ncol = K)
  #initializing parameters
  temp_par <- init_pars_latent(deg1, deg2, blocks)
  #browser()
  for (i in 1:K) {
    for (j in 1:K) {
      #browser()
      if(i > j) next # call for upper triangular part
      flag <- T
      temp_tol <- 10^(-12)
      while (flag) {
        flag <- F
        temp_sol <- tryCatch(
          expr = {
            #message(paste("Trying ", (log10(temp_tol) + 13), "th time."))
            withTimeout({
              if(is.null(lam)){
                c(estim_fn(c(i,j), A1[blocks ==i, blocks ==j], A2[blocks ==i, blocks ==j],
                                                   deg_prod1[blocks ==i, blocks ==j],
                                                   deg_prod2[blocks ==i, blocks ==j],
                                                   init_par = c(temp_par[[1]][i,j], temp_par[[2]][i,j], temp_par[[3]][i,j]),
                                                   tol_val = temp_tol, which(fold_mat[blocks ==i, blocks ==j]))[[1]])
                } else{
                  c(estim_fn(c(i,j), A1[blocks ==i, blocks ==j], A2[blocks ==i, blocks ==j],
                             deg_prod1[blocks ==i, blocks ==j],
                             deg_prod2[blocks ==i, blocks ==j],
                             init_par = c(temp_par[[1]][i,j], temp_par[[2]][i,j], temp_par[[3]][i,j]),
                             tol_val = temp_tol, which(fold_mat[blocks ==i, blocks ==j]), lam)[[1]])
                }
              },
                        timeout = max(180, 300 * (n/1000)^2 * (4/K^2)))
          },
          TimeoutException = function(ex) {
            message("Timeout Error. Retrying with new initial parameters.")
            flag <<- T
            temp_par <<- init_pars_latent(deg1, deg2, blocks)
          },
          error = function(e){
            message(e)
            message(paste("Error. Retrying ", (log10(temp_tol) + 13), "th time."))
            flag <<- T
            temp_tol <<- temp_tol * 10
            if (temp_tol > 1) {
              flag <<- F
            }
            temp_par <<- init_pars_latent(deg1, deg2, blocks)
            return(c(10, 10, 10))
          },
          warning = function(w){
            message("Warning.  Retrying. ")
            flag <<- T
            temp_par <<- init_pars_latent(deg1, deg2, blocks)
            return(c(10, 10, 10))
          },
          finally = {}
        )
      }
      sol_p1[i,j] <- temp_sol[1]
      sol_p2[i,j] <- temp_sol[2]
      sol_rho[i,j] <- temp_sol[3]
      sol_p1[j,i] <- temp_sol[1]
      sol_p2[j,i] <- temp_sol[2]
      sol_rho[j,i] <- temp_sol[3]
      #print(paste(i,j, sep = " "))
    }
  }
  return(list(sol_p1, sol_p2, sol_rho))
}

block_pair_estim_latent <- function(block_pair, A1, A2, deg_prod1, deg_prod2, 
                                    init_par = NULL, tol_val = NULL, fold = NULL){
  n <- dim(A1)[1]
  if(is.null(fold)){
    #fold stores the edge indices to be used for estimating parameters
    fold <- matrix(1, nrow = n, ncol = n)
    fold <- which(fold > 0)
  }
  #browser()
  #define objective function
  obj_fn_latent <- function(p){
    obj_sum <- 0
    temp <- 0
    #browser()
    #solnl first computes objective function and then checks constraints. 
    #So some of the arguments may be negative (but v. small)
    #A special case is introduced for each case below to handle this.
    for (e in fold) {
      i <- edge_idx(e, n)[1]
      j <- edge_idx(e, n)[2]
      if(block_pair[1] == block_pair[2] & i >= j) next #consider only upper triangular part
      if (p[3] >= 0) {
        temp <- p[3] * (min(p[1] * deg_prod1[i,j], p[2] * deg_prod2[i,j]) - 
                          p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j]) 
      } else{
        temp <- abs(p[3]) * (max(p[1] * deg_prod1[i,j] + p[2] * deg_prod2[i,j] - 1, 0) - 
                          p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j]) 
      }
      if(A1[i,j] == 1 & A2[i,j] == 0){
        if(p[1] * deg_prod1[i,j] * (1 - p[2] * deg_prod2[i,j]) - temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log(p[1] * deg_prod1[i,j] * (1 - p[2] * deg_prod2[i,j]) - temp)
      }
      if(A1[i,j] == 0 & A2[i,j] == 1){
        if(p[2] * deg_prod2[i,j] * (1 - p[1] * deg_prod1[i,j]) - temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log(p[2] * deg_prod2[i,j] * (1 - p[1] * deg_prod1[i,j]) - temp)
      }
      if(A1[i,j] == 1 & A2[i,j] == 1){
        if(p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j] + temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum = obj_sum + log(p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j] + temp)
      }
      if(A1[i,j] == 0 & A2[i,j] == 0){
        if((1 - p[1] * deg_prod1[i,j]) * (1 - p[2] * deg_prod2[i,j]) + temp  <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log((1 - p[1] * deg_prod1[i,j]) * (1 - p[2] * deg_prod2[i,j]) + temp)
      }
    }
    #browser()
    return(-1 * obj_sum)
  }
  #define linear inequality constraints
  A_ineq <- NULL
  B_ineq <- NULL
  #p1, p2 >=0, rho >= -1
  A_ineq <- rbind(A_ineq, diag(1, nrow = 3, ncol = 3))
  B_ineq <- c(B_ineq, c(0,0,-1))
  #upper bounds on p1, p2, q
  A_ineq <- rbind(A_ineq, c(-1, 0, 0))
  B_ineq <- c(B_ineq, -1/max(deg_prod1))
  A_ineq <- rbind(A_ineq, c(0, -1, 0))
  B_ineq <- c(B_ineq, -1/max(deg_prod2))
  A_ineq <- rbind(A_ineq, c(0, 0, -1))
  B_ineq <- c(B_ineq, -1)
  A_ineq <- -1 * A_ineq
  B_ineq <- -1 * B_ineq
  #find optimum parameters
  if(is.null(init_par)){
    init_par <- c(0.2, 0.2, 0.04)
  }
  if(is.null(tol_val)) tol_val <- 10^(-12)
  opt_par <- solnl(X = init_par, objfun = obj_fn_latent, A = A_ineq, B = B_ineq,
                   tolFun = tol_val, tolX = tol_val, maxnFun = 1000)
  return(opt_par)
}

block_pair_estim_lat_ot <- function(block_pair, A1, A2, deg_prod1, deg_prod2, 
                                    init_par = NULL, tol_val = NULL, fold = NULL, lam = 0.5){
  n <- dim(A1)[1]
  if(is.null(fold)){
    #fold stores the edge indices to be used for estimating parameters
    fold <- matrix(1, nrow = n, ncol = n)
    fold <- which(fold > 0)
  }
  #browser()
  #define objective function
  obj_fn_latent <- function(p){
    obj_sum <- 0
    temp <- 0
    #browser()
    #solnl first computes objective function and then checks constraints. 
    #So some of the arguments may be negative (but v. small)
    #A special case is introduced for each case below to handle this.
    for (e in fold) {
      i <- edge_idx(e, n)[1]
      j <- edge_idx(e, n)[2]
      if(block_pair[1] == block_pair[2] & i >= j) next #consider only upper triangular part
      if (p[3] >= 0) {
        temp <- p[3] * (min(p[1] * deg_prod1[i,j], p[2] * deg_prod2[i,j]) - 
                          p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j]) 
      } else{
        temp <- abs(p[3]) * (max(p[1] * deg_prod1[i,j] + p[2] * deg_prod2[i,j] - 1, 0) - 
                               p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j]) 
      }
      if(A1[i,j] == 1 & A2[i,j] == 0){
        if(p[1] * deg_prod1[i,j] * (1 - p[2] * deg_prod2[i,j]) - temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log(p[1] * deg_prod1[i,j] * (1 - p[2] * deg_prod2[i,j]) - temp)
      }
      if(A1[i,j] == 0 & A2[i,j] == 1){
        if(p[2] * deg_prod2[i,j] * (1 - p[1] * deg_prod1[i,j]) - temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log(p[2] * deg_prod2[i,j] * (1 - p[1] * deg_prod1[i,j]) - temp)
      }
      if(A1[i,j] == 1 & A2[i,j] == 1){
        if(p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j] + temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum = obj_sum + log(p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j] + temp)
      }
      if(A1[i,j] == 0 & A2[i,j] == 0){
        if((1 - p[1] * deg_prod1[i,j]) * (1 - p[2] * deg_prod2[i,j]) + temp  <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log((1 - p[1] * deg_prod1[i,j]) * (1 - p[2] * deg_prod2[i,j]) + temp)
      }
      obj_sum <- obj_sum - lam * abs(p[1] * deg_prod1[i,j] -  p[2] * deg_prod2[i,j])
    }
    #browser()
    return(-1 * obj_sum)
  }
  #define linear inequality constraints
  A_ineq <- NULL
  B_ineq <- NULL
  #p1, p2 >=0, rho >= -1
  A_ineq <- rbind(A_ineq, diag(1, nrow = 3, ncol = 3))
  B_ineq <- c(B_ineq, c(0,0,-1))
  #upper bounds on p1, p2, q
  A_ineq <- rbind(A_ineq, c(-1, 0, 0))
  B_ineq <- c(B_ineq, -1/max(deg_prod1))
  A_ineq <- rbind(A_ineq, c(0, -1, 0))
  B_ineq <- c(B_ineq, -1/max(deg_prod2))
  A_ineq <- rbind(A_ineq, c(0, 0, -1))
  B_ineq <- c(B_ineq, -1)
  A_ineq <- -1 * A_ineq
  B_ineq <- -1 * B_ineq
  #find optimum parameters
  if(is.null(init_par)){
    init_par <- c(0.2, 0.2, 0.04)
  }
  if(is.null(tol_val)) tol_val <- 10^(-12)
  opt_par <- solnl(X = init_par, objfun = obj_fn_latent, A = A_ineq, B = B_ineq,
                   tolFun = tol_val, tolX = tol_val, maxnFun = 1000)
  return(opt_par)
}

block_pair_estim_lat_pen <- function(block_pair, A1, A2, deg_prod1, deg_prod2, 
                                    init_par = NULL, tol_val = NULL, fold = NULL, lam = 0.5){
  n <- dim(A1)[1]
  if(is.null(fold)){
    #fold stores the edge indices to be used for estimating parameters
    fold <- matrix(1, nrow = n, ncol = n)
    fold <- which(fold > 0)
  }
  #browser()
  #define objective function
  obj_fn_latent <- function(p){
    obj_sum <- 0
    temp <- 0
    #browser()
    #solnl first computes objective function and then checks constraints. 
    #So some of the arguments may be negative (but v. small)
    #A special case is introduced for each case below to handle this.
    for (e in fold) {
      i <- edge_idx(e, n)[1]
      j <- edge_idx(e, n)[2]
      if(block_pair[1] == block_pair[2] & i >= j) next #consider only upper triangular part
      if (p[3] >= 0) {
        temp <- p[3] * (min(p[1] * deg_prod1[i,j], p[2] * deg_prod2[i,j]) - 
                          p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j]) 
      } else{
        temp <- abs(p[3]) * (max(p[1] * deg_prod1[i,j] + p[2] * deg_prod2[i,j] - 1, 0) - 
                               p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j]) 
      }
      if(A1[i,j] == 1 & A2[i,j] == 0){
        if(p[1] * deg_prod1[i,j] * (1 - p[2] * deg_prod2[i,j]) - temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log(p[1] * deg_prod1[i,j] * (1 - p[2] * deg_prod2[i,j]) - temp)
      }
      if(A1[i,j] == 0 & A2[i,j] == 1){
        if(p[2] * deg_prod2[i,j] * (1 - p[1] * deg_prod1[i,j]) - temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log(p[2] * deg_prod2[i,j] * (1 - p[1] * deg_prod1[i,j]) - temp)
      }
      if(A1[i,j] == 1 & A2[i,j] == 1){
        if(p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j] + temp <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum = obj_sum + log(p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j] + temp)
      }
      if(A1[i,j] == 0 & A2[i,j] == 0){
        if((1 - p[1] * deg_prod1[i,j]) * (1 - p[2] * deg_prod2[i,j]) + temp  <= 0){
          obj_sum <- -Inf
          break
        }
        obj_sum <- obj_sum + log((1 - p[1] * deg_prod1[i,j]) * (1 - p[2] * deg_prod2[i,j]) + temp)
      }
      obj_sum <- obj_sum - lam * (p[1] * deg_prod1[i,j] +  p[2] * deg_prod2[i,j] - 
                                    2*(p[1] * deg_prod1[i,j] * p[2] * deg_prod2[i,j]) -
                                    2 * temp)
    }
    #browser()
    return(-1 * obj_sum)
  }
  #define linear inequality constraints
  A_ineq <- NULL
  B_ineq <- NULL
  #p1, p2 >=0, rho >= -1
  A_ineq <- rbind(A_ineq, diag(1, nrow = 3, ncol = 3))
  B_ineq <- c(B_ineq, c(0,0,-1))
  #upper bounds on p1, p2, q
  A_ineq <- rbind(A_ineq, c(-1, 0, 0))
  B_ineq <- c(B_ineq, -1/max(deg_prod1))
  A_ineq <- rbind(A_ineq, c(0, -1, 0))
  B_ineq <- c(B_ineq, -1/max(deg_prod2))
  A_ineq <- rbind(A_ineq, c(0, 0, -1))
  B_ineq <- c(B_ineq, -1)
  A_ineq <- -1 * A_ineq
  B_ineq <- -1 * B_ineq
  #find optimum parameters
  if(is.null(init_par)){
    init_par <- c(0.2, 0.2, 0.04)
  }
  if(is.null(tol_val)) tol_val <- 10^(-12)
  opt_par <- solnl(X = init_par, objfun = obj_fn_latent, A = A_ineq, B = B_ineq,
                   tolFun = tol_val, tolX = tol_val, maxnFun = 1000)
  return(opt_par)
}

gen_auc_corr_dcsbm <- function(A1, A2, deg1, deg2, blocks, dcsbm_estim = block_pair_estim_latent, lam = NULL){
  folds <- createFolds(which(upper.tri(A1)), 5, list = T, returnTrain = F) #for cross-validation
  folds <- lapply(folds, function(x) which(upper.tri(A1))[x])
  #browser()
  auc_sum <- 0
  roc_list <- list()
  pars_list <- list()
  for (i in 1:5) {
    #print(i)
    #browser()
    if(is.null(lam)){
      pars_list[[i]] <- estim_dcsbm_latent(A1, A2, deg1, deg2, blocks, setdiff(which(upper.tri(A1)), folds[[i]]), estim = dcsbm_estim)
    } else{
      pars_list[[i]] <- estim_dcsbm_latent(A1, A2, deg1, deg2, blocks, setdiff(which(upper.tri(A1)), folds[[i]]), estim = dcsbm_estim, lam)
    }
    roc_list[[i]] <- compute_tpr_fpr_dcsbm(A1, A2, deg1, deg2, blocks, pars_list[[i]], folds[[i]])
  }
  for (i in 1:5) {    auc_sum <- auc_sum + auc(roc_list[[i]])  }
  return(list(pars_list, roc_list, auc_sum/5))
}

compute_tpr_fpr_dcsbm <- function(A1, A2, deg1, deg2, blocks, pars_list, fold){
  n <- dim(A1)[1]
  e_id <- c(0,0)
  block_id <- c(0,0)
  temp_par <- rep(0,4)
  A1_e <- 0
  pred_pars <- rep(0, length(fold))
  counter <- 1
  for (e in fold) {
    e_id <- edge_idx(e, n)
    assert_that(e_id[1] < e_id[2])
    if(e_id[1] >= e_id[2]) next #remove this line later. Take note of var counter
    block_id <- blocks[e_id]
    temp_par <- sapply(pars_list, function(x){ 
      x[block_id[1], block_id[2]]
    })
    temp_par[1] <- temp_par[1] * deg1[e_id[1]] * deg1[e_id[2]]
    temp_par[2] <- temp_par[2] * deg2[e_id[1]] * deg2[e_id[2]]
    if(temp_par[3] >= 0){
      temp_par[3] <- temp_par[3] * (min(temp_par[1], temp_par[2]) -
                                      (temp_par[1] * temp_par[2]))
    } else{
      temp_par[3] <- abs(temp_par[3]) * (max(temp_par[1] + temp_par[2] - 1, 0) -
                                           (temp_par[1] * temp_par[2]))
    }
    A1_e <- A1[e_id[1], e_id[2]]
    if(A1_e == 0){
        pred_pars[counter] <- (temp_par[2] * (1 - temp_par[1]) - temp_par[3])/
                              (1 - temp_par[1]) #(p2-q)/(1-p1)
    } else{
      if(temp_par[1] == 0) pred_pars[counter] <- 0
      else{
        pred_pars[counter] <- (temp_par[1] * temp_par[2] + temp_par[3])/
                            (temp_par[1]) #q/p1
      }
    }
    counter <- counter + 1
  }
  #browser()
  #assert_that(sum(is.na(pred_pars)) == 0)
  if(sum(is.na(pred_pars)) > 0) browser()
  assert_that(counter == (length(fold) + 1))
  roc_dcsbm <- roc(A2[fold], pred_pars)
  return(roc_dcsbm)
}


#some code for testing out the functions
symm_mat <- function(X, normalize = NULL){
  X[lower.tri(X, diag = F)] <- 0
  X <- (X + t(X))
  diag(X) <- diag(X)/2
  return(X)
}

node_gp <- function(i, n, K){
  s <- floor(n/K) # s is block size
  return(min((i-1)%/%s + 1, K)) 
}

###########################
#code for testing out
n <- 100
K <- 2
dir_var <- 1
p1 <- matrix(0, nrow = K, ncol = K)
#p1 <- matrix(0.03, nrow = K, ncol = K)
p1[upper.tri(p1, diag = T)] <- c(rdirichlet(1, rep(dir_var, K*(K+1)/2)))/n * K*(K+1)/2 
p2 <- matrix(0, nrow = K, ncol = K)
p2[upper.tri(p2, diag = T)] <- c(rdirichlet(1, rep(dir_var, K*(K+1)/2)))/n * K*(K+1)/2 
#p2 <- matrix(0.04, nrow = K, ncol = K)
p1 <- symm_mat(p1)
p2 <- symm_mat(p2)
s_deg1 <- pow_law_sample(n, 1, 10, 2)
s_deg2 <- pow_law_sample(n, 1, 10, 2)
#s_deg1 <- c(rep(1, n/2), rep(2, n/2))
#s_deg2 <- c(rep(1, n/2), rep(2, n/2))
#s_deg1 <- rep(1, n)
#s_deg2 <- rep(1, n)
p2 <- p1
s_deg2 <- s_deg1
temp_blocks <- sapply(1:n, function(x) node_gp(x, n, K))
temp_rho <- matrix(-0.95, nrow = K, ncol = K) 
assert_that(sum(abs(temp_rho - t(temp_rho))) == 0)
check_feas_par(p1, p2, temp_rho, s_deg1, s_deg2, temp_blocks)
nets <- corr_dcsbm_latent(p1, p2, temp_rho, s_deg1, s_deg2, temp_blocks)
out_par <- estim_dcsbm_latent(nets[[1]], nets[[2]], s_deg1, s_deg2, temp_blocks)
#print(c(out_par[[1]], out_par[[2]], out_par[[3]]))
print(out_par)
out_auc <- gen_auc_corr_dcsbm(nets[[1]], nets[[2]], s_deg1, s_deg2, temp_blocks)
out_auc[[3]]
out_auc <- gen_auc_corr_dcsbm(nets[[1]], nets[[2]], s_deg1, s_deg2, temp_blocks, dcsbm_estim = block_pair_estim_lat_pen, lam = 1)
out_auc[[3]]


#write pars to file
temp_dat <- list(nets, temp_blocks)
saveRDS(temp_dat, file = "temp_output.rds")
#ml pars in indep case
fld <- NULL
my_edges <- get_edges_by_blocks_sbm(nets[[1]], nets[[2]], temp_blocks, fld)
ml_p1 <- (my_edges[[1]] + my_edges[[3]])/my_edges$total
ml_p2 <- (my_edges[[2]] + my_edges[[3]])/my_edges$total
ml_p3 <- (my_edges[[3]])/my_edges$total


#####################################################
#Code for output
n <- 100
K <- 2
dir_var <- 1
p1 <- matrix(0, nrow = K, ncol = K)
p1[upper.tri(p1, diag = T)] <- c(rdirichlet(1, rep(dir_var, K*(K+1)/2)))/n * K*(K+1)
p2 <- matrix(0, nrow = K, ncol = K)
p2[upper.tri(p2, diag = T)] <- c(rdirichlet(1, rep(dir_var, K*(K+1)/2)))/n * K*(K+1) 
p1 <- symm_mat(p1)
p2 <- symm_mat(p2)
s_deg1 <- pow_law_sample(n, 1, 10, 2)
s_deg2 <- pow_law_sample(n, 1, 10, 2)
#p1 <- matrix(0.025, nrow = 1, ncol = 1)
#s_deg2 <- s_deg1
temp_blocks <- sapply(1:n, function(x) node_gp(x, n, K))
t_list <- 0:10
out_auc <- matrix(0, nrow = 2, ncol = 11)
for (t in 4:10) {
  #pos. correlated
  for (j in 1:5) {
    temp_rho <- matrix(t/10, nrow = K, ncol = K)
    check_feas_par(p1, p2, temp_rho, s_deg1, s_deg2, temp_blocks)
    nets <- corr_dcsbm_latent(p1, p2, temp_rho, s_deg1, s_deg2, temp_blocks)
    saveRDS(file = paste0(t,j,"pos"), object = list(p1, p2, temp_rho, nets))
    temp_auc <- gen_auc_corr_dcsbm(nets[[1]], nets[[2]], s_deg1, s_deg2, temp_blocks)
    out_auc[1, t+1] <- out_auc[1, t+1] + temp_auc[[3]]
  }
  out_auc[1, t+1] <- out_auc[1, t+1]/5
  #neg. correlated
  for (j in 1:5) {
    temp_rho <- matrix(-t/10, nrow = K, ncol = K)
    check_feas_par(p1, p2, temp_rho, s_deg1, s_deg2, temp_blocks)
    nets <- corr_dcsbm_latent(p1, p2, temp_rho, s_deg1, s_deg2, temp_blocks)
    saveRDS(file = paste0(t,j,"neg"), object = list(p1, p2, temp_rho, nets))
    temp_auc <- gen_auc_corr_dcsbm(nets[[1]], nets[[2]], s_deg1, s_deg2, temp_blocks)
    out_auc[2, t+1] <- out_auc[2, t+1] + temp_auc[[3]]
  }
  out_auc[2, t+1] <- out_auc[2, t+1]/5
}
saveRDS(out_auc, file = "AUC_list.rds")
#with ot penalty, lam=0.5
t_list <- 0:10
out_auc <- matrix(0, nrow = 2, ncol = 11)
for (t in t_list) {
  #pos. correlated
  for (j in 1:5) {
    tpars <- readRDS(file = paste0(t, j, "pos"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_ot,
                                   lam = 0.5)
    out_auc[1, t+1] <- out_auc[1, t+1] + temp_auc[[3]]
  }
  out_auc[1, t+1] <- out_auc[1, t+1]/5
  #neg. correlated
  for (j in 1:5) {
    tpars <- readRDS(file = paste0(t, j, "neg"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_ot,
                                   lam = 0.5)
    out_auc[2, t+1] <- out_auc[2, t+1] + temp_auc[[3]]
  }
  out_auc[2, t+1] <- out_auc[2, t+1]/5
}
saveRDS(out_auc, file = "AUC_list_OT_5.rds")

#with ot penalty, lam=1
t_list <- 0:10
out_auc <- matrix(0, nrow = 2, ncol = 11)
for (t in 8:10) {
  #pos. correlated
  for (j in 1:5) {
    print(paste(t, j, "pos"))
    tpars <- readRDS(file = paste0(t, j, "pos"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_ot,
                                   lam = 1)
    print(paste(t, j, "pos"))
    out_auc[1, t+1] <- out_auc[1, t+1] + temp_auc[[3]]
  }
  out_auc[1, t+1] <- out_auc[1, t+1]/5
  #neg. correlated
  for (j in 1:5) {
    print(paste(t, j, "neg"))
    tpars <- readRDS(file = paste0(t, j, "neg"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_ot,
                                   lam = 1)
    out_auc[2, t+1] <- out_auc[2, t+1] + temp_auc[[3]]
  }
  out_auc[2, t+1] <- out_auc[2, t+1]/5
}
saveRDS(out_auc, file = "AUC_list_OT_10.rds")

###############
#code with mass tran. penalty
# lam=0.5
t_list <- 0:10
out_auc <- matrix(0, nrow = 2, ncol = 11)
for (t in 5:10) {
  #pos. correlated
  for (j in 1:5) {
    print(paste(t, j, "pos"))
    tpars <- readRDS(file = paste0(t, j, "pos"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_pen,
                                   lam = 0.5)
    out_auc[1, t+1] <- out_auc[1, t+1] + temp_auc[[3]]
  }
  out_auc[1, t+1] <- out_auc[1, t+1]/5
  #neg. correlated
  for (j in 1:5) {
    print(paste(t, j, "neg"))
    tpars <- readRDS(file = paste0(t, j, "neg"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_pen,
                                   lam = 0.5)
    out_auc[2, t+1] <- out_auc[2, t+1] + temp_auc[[3]]
  }
  out_auc[2, t+1] <- out_auc[2, t+1]/5
}
saveRDS(out_auc, file = "AUC_list_PEN_5.rds")

#with ot penalty, lam=2
t_list <- 0:10
out_auc <- matrix(0, nrow = 2, ncol = 11)
for (t in 8:10) {
  #pos. correlated
  for (j in 1:5) {
    print(paste(t, j, "pos"))
    tpars <- readRDS(file = paste0(t, j, "pos"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_pen,
                                   lam = 2)
    out_auc[1, t+1] <- out_auc[1, t+1] + temp_auc[[3]]
  }
  out_auc[1, t+1] <- out_auc[1, t+1]/5
  #neg. correlated
  for (j in 1:5) {
    print(paste(t, j, "neg"))
    tpars <- readRDS(file = paste0(t, j, "neg"))
    check_feas_par(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    nets <- corr_dcsbm_latent(tpars[[1]], tpars[[2]], tpars[[3]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks)
    temp_auc <- gen_auc_corr_dcsbm(tpars[[4]][[1]], tpars[[4]][[2]], tpars[[4]][[3]], tpars[[4]][[4]], temp_blocks, 
                                   dcsbm_estim = block_pair_estim_lat_pen,
                                   lam = 2)
    out_auc[2, t+1] <- out_auc[2, t+1] + temp_auc[[3]]
  }
  out_auc[2, t+1] <- out_auc[2, t+1]/5
}
saveRDS(out_auc, file = "AUC_list_PEN_20.rds")

#########################################
#Plot to view PEN output
out_auc <- readRDS(file = "AUC_list.rds")
temp_df1 <- matrix(0, nrow = 21, ncol = 2)
temp_df1[,1] <- c(-1* rev(1:10)/10, (0:10)/10)
temp_df1[,2] <- c(rev(out_auc[2, 2:11]), out_auc[1, 1:11])
temp_df1 <- as.data.frame(temp_df1)
colnames(temp_df1) <- c("t", "AUC")

out_auc <- readRDS(file = "AUC_list_PEN_5.rds")
temp_df2 <- matrix(0, nrow = 21, ncol = 2)
temp_df2[,1] <- c(-1* rev(1:10)/10, (0:10)/10)
temp_df2[,2] <- c(rev(out_auc[2, 2:11]), out_auc[1, 1:11])
temp_df2 <- as.data.frame(temp_df2)
colnames(temp_df2) <- c("t", "AUC")

out_auc <- readRDS(file = "AUC_list_PEN_20.rds")
temp_df3 <- matrix(0, nrow = 21, ncol = 2)
temp_df3[,1] <- c(-1* rev(1:10)/10, (0:10)/10)
temp_df3[,2] <- c(rev(out_auc[2, 2:11]), out_auc[1, 1:11])
temp_df3 <- as.data.frame(temp_df3)
colnames(temp_df3) <- c("t", "AUC")

temp_df1$lam <- rep(0, dim(temp_df1)[1])
temp_df2$lam <- rep(0.5, dim(temp_df2)[1])
temp_df3$lam <- rep(2, dim(temp_df3)[1])

temp_df <- rbind(temp_df1, temp_df2, temp_df3)
temp_df$lam <- as.factor(temp_df$lam)
ggplot(temp_df, aes(t, AUC)) + geom_line(aes(color = lam))

#################################
#Plot to view OT output
out_auc <- readRDS(file = "AUC_list.rds")
temp_df1 <- matrix(0, nrow = 21, ncol = 2)
temp_df1[,1] <- c(-1* rev(1:10)/10, (0:10)/10)
temp_df1[,2] <- c(rev(out_auc[2, 2:11]), out_auc[1, 1:11])
temp_df1 <- as.data.frame(temp_df1)
colnames(temp_df1) <- c("t", "AUC")

out_auc <- readRDS(file = "AUC_list_OT_5.rds")
temp_df2 <- matrix(0, nrow = 21, ncol = 2)
temp_df2[,1] <- c(-1* rev(1:10)/10, (0:10)/10)
temp_df2[,2] <- c(rev(out_auc[2, 2:11]), out_auc[1, 1:11])
temp_df2 <- as.data.frame(temp_df2)
colnames(temp_df2) <- c("t", "AUC")

out_auc <- readRDS(file = "AUC_list_OT_10.rds")
temp_df3 <- matrix(0, nrow = 21, ncol = 2)
temp_df3[,1] <- c(-1* rev(1:10)/10, (0:10)/10)
temp_df3[,2] <- c(rev(out_auc[2, 2:11]), out_auc[1, 1:11])
temp_df3 <- as.data.frame(temp_df3)
colnames(temp_df3) <- c("t", "AUC")

temp_df1$lam <- rep(0, dim(temp_df1)[1])
temp_df2$lam <- rep(0.5, dim(temp_df2)[1])
temp_df3$lam <- rep(1, dim(temp_df3)[1])

temp_df <- rbind(temp_df1, temp_df2, temp_df3)
temp_df$lam <- as.factor(temp_df$lam)
ggplot(temp_df, aes(t, AUC)) + geom_line(aes(color = lam))


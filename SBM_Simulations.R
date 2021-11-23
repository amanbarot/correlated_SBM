#Simulations for correlated SBMs
library(CVXR)
library(ggplot2)
library(reshape2)
library(cowplot)

#Compare blocks models with 3 blocks.
#Define model parameters
K <- 3
real_p1 <- matrix(c(0.3,0.15,0.15, 0.15,0.4,0.15, 0.15,0.15,0.5), nrow = K, byrow = T)
real_p2 <- matrix(c(0.3,0.15,0.15, 0.15,0.3,0.15, 0.15,0.15,0.4), nrow = K, byrow = T)
real_q <- matrix(0.09, nrow = K, ncol = K)

real_p10 <- real_p1 - real_q
real_p01 <- real_p2 - real_q
real_p11 <- real_q
real_p00 <- 1 - real_p1 - real_p2 + real_q

n <- 99
N <- n*(n-1)/2
block_n <- n/K 
block_N <- block_n * (block_n - 1)/2

my_edges <- array(dim = c(K,K,4))

#testing pars
ntests <- 1 #tests per lambda
nlam <- c(((1:5) - 1)/10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10) 

#Solve with OT regularization
ot_out <- array(0, dim = c(length(nlam), K, K, 4))

for (i in 1:ntests) {
  #Generate edges
  for (k in 1:K) {
    for (l in 1:K) {
      if(k <=l){
      my_edges[k,l,] <- rmultinom(1, block_N, c(real_p10[k,l], real_p01[k,l], real_p11[k,l],
                                                real_p00[k,l]))
      }
      else{
        my_edges[k,l,] <- my_edges[l,k,]
      }
    }
  }
  #Vary lambda, the reg parameter
  for (j in 1:length(nlam)) {
    lam <- nlam[j]
    p1 <- Variable(K,K)
    p2 <- Variable(K,K)
    p3 <- Variable(1)
    p4 <- Variable(K,K)
    objective <- Minimize(sum(-log(p1) * my_edges[,,1] - log(p2) * my_edges[,,2] -
                                log(p3) * my_edges[,,3] - log(p4) * my_edges[,,4] +
                                lam * block_N * (p1 + p2 - 2*min(p1, p2))))
    my_constraints <- list(p1>=0, p1<=1, p2>=0, p2<=1, p3>=0, p3<=1, p4>=0, p4<=1,
                           p1+p2+p3+p4==1, p1-t(p1)==0, p2-t(p2)==0, p4-t(p4)==0)
    my_prob <- Problem(objective, my_constraints)
    sol_ot <- solve(my_prob)
    ot_out[j,,,1] <- ot_out[j,,,1] + sol_ot$getValue(p1)
    ot_out[j,,,2] <- ot_out[j,,,2] + sol_ot$getValue(p2)
    ot_out[j,,,3] <- ot_out[j,,,3] + sol_ot$getValue(p3)
    ot_out[j,,,4] <- ot_out[j,,,4] + sol_ot$getValue(p4)
  }
}
ot_out <- ot_out/ntests

#Visualize p1,p2,p3,p4
temp_df <- as.data.frame(matrix(0, nrow = length(nlam) * K*K * 4, ncol = 4))
colnames(temp_df) <- c("lambda", "block_idx", "par", "value")
temp_df$block_idx <- as.character(temp_df$block_idx)
temp_df$par <- as.character(temp_df$par)
#temp_df$value <- as.character(temp_df$value)
counter <- 1
for (i in 1:length(nlam)) {
  for (j in 1:(K*K)) {
    for (k in 1:4) {
      temp_df[counter,1] <- nlam[i]
      temp_df[counter,2] <- paste((j-1)%/%K + 1, (j-1)%%K + 1, sep = ",")
      if(k==1){
        temp_df[counter,3] <- paste("p10")
      }
      if(k==2){
        temp_df[counter,3] <- paste("p01")
      }
      if(k==3){
        temp_df[counter,3] <- paste("p11")
      }
      if(k==4){
        temp_df[counter,3] <- paste("p00")
      }
      temp_df[counter,4] <- ot_out[i, (j-1)%/%K + 1, (j-1)%%K + 1, k]
      #print(paste(i,j,k))
      #print(temp_df[counter,])
      counter <- counter + 1
    }
  }
}
ggplot(data = temp_df, aes(lambda, value, group = par, color = par)) +
  geom_point() + facet_wrap(~block_idx, ncol = 3)

#################################################
#Solve by penalization based on the candidate coupling
coup_out <- array(0, dim = c(length(nlam), K, K, 4))

for (i in 1:ntests) {
  #Generate edges
  for (k in 1:K) {
    for (l in 1:K) {
      if(k <=l){
        my_edges[k,l,] <- rmultinom(1, block_N, c(real_p10[k,l], real_p01[k,l], real_p11[k,l],
                                                  real_p00[k,l]))
      }
      else{
        my_edges[k,l,] <- my_edges[l,k,]
      }
    }
  }
  #Vary lambda, the reg parameter
  for (j in 1:length(nlam)) {
    lam <- nlam[j]
    p1 <- Variable(K,K)
    p2 <- Variable(K,K)
    p3 <- Variable(1)
    p4 <- Variable(K,K)
    objective <- Minimize(sum(-log(p1) * my_edges[,,1] - log(p2) * my_edges[,,2] -
                                log(p3) * my_edges[,,3] - log(p4) * my_edges[,,4] +
                                lam * block_N * (p1 + p2)))
    my_constraints <- list(p1>=0, p1<=1, p2>=0, p2<=1, p3>=0, p3<=1, p4>=0, p4<=1,
                           p1+p2+p3+p4==1, p1-t(p1)==0, p2-t(p2)==0, p4-t(p4)==0)
    my_prob <- Problem(objective, my_constraints)
    sol_ot <- solve(my_prob)
    coup_out[j,,,1] <- coup_out[j,,,1] + sol_ot$getValue(p1)
    coup_out[j,,,2] <- coup_out[j,,,2] + sol_ot$getValue(p2)
    coup_out[j,,,3] <- coup_out[j,,,3] + sol_ot$getValue(p3)
    coup_out[j,,,4] <- coup_out[j,,,4] + sol_ot$getValue(p4)
  }
}
coup_out <- coup_out/ntests

#Visualize p1,p2,p3,p4
temp_df <- as.data.frame(matrix(0, nrow = length(nlam) * K*K * 4, ncol = 4))
colnames(temp_df) <- c("lambda", "block_idx", "par", "value")
temp_df$block_idx <- as.character(temp_df$block_idx)
temp_df$par <- as.character(temp_df$par)
counter <- 1
for (i in 1:length(nlam)) {
  for (j in 1:(K*K)) {
    for (k in 1:4) {
      temp_df[counter,1] <- nlam[i]
      temp_df[counter,2] <- paste((j-1)%/%K + 1, (j-1)%%K + 1, sep = ",")
      if(k==1){
        temp_df[counter,3] <- paste("p10")
      }
      if(k==2){
        temp_df[counter,3] <- paste("p01")
      }
      if(k==3){
        temp_df[counter,3] <- paste("p11")
      }
      if(k==4){
        temp_df[counter,3] <- paste("p00")
      }
      temp_df[counter,4] <- coup_out[i, (j-1)%/%K + 1, (j-1)%%K + 1, k]
      counter <- counter + 1
    }
  }
}
ggplot(data = temp_df, aes(lambda, value, group = par, color = par)) +
  geom_point() + facet_wrap(~block_idx, ncol = 3)

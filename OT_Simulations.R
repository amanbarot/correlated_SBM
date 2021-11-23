#Simulations for correlated Erdos-Renyi graphs
library(CVXR)
library(ggplot2)
library(reshape2)
library(cowplot)

real_p <- c(0.3, 0.3)
real_q <- 0.1
n <- 100
N <- n*(n-1)/2
n_lam <- 21 #11 
ntests <- 10 #10
out_p_ot <- matrix(0, nrow = n_lam, ncol = 4)
#out_q_ot <- rep(0, n_lam)
out_p_pen <- matrix(0, nrow = n_lam, ncol = 4)
#out_q_pen <- rep(0, n_lam)
for (i in 1:ntests) {
  my_edges <- c(rmultinom(1, N, c(real_p[1] - real_q, real_p[2] - real_q, 
                                  real_q, 1 - real_p[1] - real_p[2] + real_q))) #p1-q, p2-q, q, 1-p1-p2+q
  p <- Variable(4)#Var for cvx solver
  #q <- Variable(1)
  for (j in 1:n_lam) {
    lam <- (j-1)/10
    #Minimize -log L + lam * OT
    objective <- Minimize(- t(my_edges) %*% log(p)
                          + lam * N * abs(p[1] - p[2])) #p1 -q -(p2-q) = p1-p2
    constraints <- list(p>=0, p<=1, sum(p)==1)
    my_prob <- Problem(objective, constraints)
    sol_ot <- solve(my_prob)
    out_p_ot[j,] <- out_p_ot[j,] + c(sol_ot$getValue(p))
    
    #Compare with penalizing with penalty for the plan
    objective <- Minimize(- t(my_edges) %*% log(p)
                          + lam * N * (p[1] + p[2])) #p1 -q +(p2-q) = p1+p2-2q
    #+ lam * N * (p[1] + p[2] - 2* q))
    constraints <- list(p>=0, p<=1, sum(p)==1)
    #constraints <- list(q >= 0, q <= 1, q <= min(p), sum(p) <= 1 + q)
    my_prob <- Problem(objective, constraints)
    sol_pen <- solve(my_prob)
    out_p_pen[j,] <- out_p_pen[j,] + c(sol_pen$getValue(p))
  }
}
#plots for ot
#plot p for 10,01,11,00
temp_df <- as.data.frame(out_p_ot/ntests)
colnames(temp_df) <- c("01", "10", "11", "00")
temp_df$lambda <- ((1:n_lam) - 1)/10
melt_df <- melt(temp_df, id.vars = "lambda")
plt_1 <- ggplot(melt_df, aes(lambda, value, group = variable, color = variable)) + geom_point() + 
  ggtitle(paste("p10 = ", real_p[1]-real_q, ", p01 = ", real_p[2]-real_q, ", p11 = ", real_q,
                ", p00 = ", 1 - real_p[1] - real_p[2] + real_q, sep = ""))
#plot for p1,p2,q
temp_df <- data.frame("p1" = temp_df[,1] + temp_df[,3], "p2" = temp_df[,2] + temp_df[,3], 
                      "q" = temp_df[,3])
temp_df$lambda <- ((1:n_lam) - 1)/10
melt_df <- melt(temp_df, id.vars = "lambda")
plt_2 <- ggplot(melt_df, aes(lambda,value, group = variable, color = variable)) + geom_point() + 
  ggtitle(paste("p1 =", real_p[1], ", p2 =", real_p[2], ", q = ", real_q, sep = ""))
plot_grid(plt_1, plt_2, nrow = 2, ncol = 1)

#plots for pen
#plot p for 10,01,11,00
temp_df <- as.data.frame(out_p_pen/ntests)
colnames(temp_df) <- c("01", "10", "11", "00")
temp_df$lambda <- ((1:n_lam) - 1)/10
melt_df <- melt(temp_df, id.vars = "lambda")
plt_1 <- ggplot(melt_df, aes(lambda, value, group = variable, color = variable)) + geom_point() + 
  ggtitle(paste("p10 = ", real_p[1]-real_q, ", p01 = ", real_p[2]-real_q, ", p11 = ", real_q,
                ", p00 = ", 1 - real_p[1] - real_p[2] + real_q, sep = ""))
#plot for p1,p2,q
temp_df <- data.frame("p1" = temp_df[,1] + temp_df[,3], "p2" = temp_df[,2] + temp_df[,3], 
                      "q" = temp_df[,3])
temp_df$lambda <- ((1:n_lam) - 1)/10
melt_df <- melt(temp_df, id.vars = "lambda")
plt_2 <- ggplot(melt_df, aes(lambda,value, group = variable, color = variable)) + geom_point() + 
  ggtitle(paste("p1 =", real_p[1], ", p2 =", real_p[2], ", q = ", real_q, sep = ""))

plot_grid(plt_1, plt_2, nrow = 2, ncol = 1)
#pick p1, p2 with varying p1+p2 and vary q.
#0,4,0,6,0.2
#0.2,0.4,0.1
#0.2,0.4,0
#0.3,0.3,0.1
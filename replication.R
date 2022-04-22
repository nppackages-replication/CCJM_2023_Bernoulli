# This is the replication file for Table 1 in Cattaneo, Chandak, Jansson and Ma (2022).
# Simulated conditional pdf estimates with coverage probabilities.

#install companion package
install.packages("lpcde")

#load companion package
library("lpcde")

# setting parameters
n = 5000
num_sims = 1000
ng = 20
CIsimul=2000
alpha=0.05
x = 0
y_grid = seq(0, 1.0, length.out=ng)


# initializing
bias = matrix(0L, nrow = ng, ncol = num_sims)
sd = aw_pw = aw_u = bias_rbc = sd_rbc = aw_pw_rbc = aw_u_rbc = bias

# Using parallel computations for simulations
#this will use all available cores on the computer
library(parallel)
numCores = detectCores()
# set the variable below if a specific number of cores are to be used
numCores

# function that runs each simulation
simtest = function(s){
  # simulating data
  data = mvtnorm::rmvnorm(5*n, sigma = matrix(c(2, -0.1, -0.1, 2), nrow=2))
  data = data[abs(data[,1])<=1.0,]
  data = data[abs(data[,2])<=1.0,]
  data = data[1:n, ]
  x_data = as.matrix(data[, 2])
  y_data = as.matrix(data[, 1])

  # bandwidth selection
  bw_star = lpcde::lpbwcde(y_data=y_data, x_data=x_data, y_grid=y_grid,
                           x=x, bw_type="mse-rot")$BW[,2]

  # estimating density
  est = lpcde::lpcde(y_data=y_data, x_data=x_data, y_grid=y_grid, x=x,
                     bw=bw_star, kernel_type = "epanechnikov")
  f_hat = est$Estimate[,3]

  # computing values of interest:
  # bias
  bias = (f_hat - true_val)

  # standard deviation
  sd = est$Estimate[,5]

  # average length of confidence interval
  aw_pw = 2*1.96*(sd)

  # confidence band computation
  corrMat = sweep(sweep(est$CovMat$CovMat, MARGIN=1, FUN="*",
                        STATS=1/est$Estimate[, "se"]),
                  MARGIN=2, FUN="*", STATS=1/est$Estimate[, "se"])
  normalSimu = try(
    MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                  Sigma=Matrix::nearPD(corrMat)$mat),
    silent=TRUE)
  z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                                FUN=function(x) {max(abs(x))}), 1 - alpha)
  print(z_val)
  aw_u = 2*z_val*sd

  # pointwise and uniform coverage
  ci_size =100*ifelse(abs(bias/sd<=1.96), 1, 0)
  cb_size = 100*ifelse(abs(bias/sd<=z_val), 1, 0)

  #rbc
  f_hat = est$Estimate[,4]

  # bias
  bias_rbc = (f_hat - true_val)

  # standard deviation
  sd_rbc = est$Estimate[,6]

  # RBC CI length
  aw_pw_rbc = 2*1.96*sd_rbc

  # CB computation
  corrMat = sweep(sweep(est$CovMat$CovMat_RBC, MARGIN=1, FUN="*",
                        STATS=1/est$Estimate[, "se_RBC"]),
                  MARGIN=2, FUN="*", STATS=1/est$Estimate[, "se_RBC"])
  normalSimu = try(
    MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                  Sigma=Matrix::nearPD(corrMat)$mat),
    silent=TRUE)
  z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                                FUN=function(x) {max(abs(x))}), 1 - alpha)
  print(z_val)
  aw_u_rbc = 2*z_val*sd_rbc

  # pointwise and uniform coverage
  ci_size_rbc = 100*ifelse((abs(bias_rbc/sd_rbc)<=1.96), 1, 0)
  cb_size_rbc = 100*ifelse((abs(bias_rbc/sd_rbc)<=z_val), 1, 0)

  # combining all results of simulation
  results = matrix(c(bw_star, bias, sd, ci_size, cb_size, aw_pw, aw_u,
                     bw_star, bias_rbc, sd_rbc, ci_size_rbc, cb_size_rbc,
                     aw_pw_rbc, aw_u_rbc),
                   ncol = 14)
  return(results)
}

pw_table = matrix(ncol=7)
for (x in c(0, 0.8, 1.0)){
  # true density function
  true_val = stats::dnorm(y_grid, mean = -0.1*x, sd=sqrt(1-0.1^2))
  # running simulations
  object = parallel::mclapply(1:num_sims, simtest, mc.cores = numCores)

  # combining results of all simulations
  avg = matrix(0L, nrow = num_sims, ncol = 14)
  for (l in 1:num_sims){
    avg[l,] = colSums(object[[l]])/ng
  }
  table_data = matrix((colSums(avg)/num_sims)[1:7], nrow=1)
  table_data = rbind(table_data, (colSums(avg)/num_sims)[8:14])

  pw_table = rbind(pw_table, table_data)
  cat("\n")
}
pw_table= pw_table[-1, ]
pw_table[, 2] = abs(pw_table[, 2])
colnames(pw_table) = c("BW", "bias", "sd", "PW Coverage", "UCB Coverage",
                         "CI AL", "CB AW")
rownames(pw_table) = c("WBC", "RBC","WBC", "RBC","WBC", "RBC")

pw_table





R = 1000 # number of MC data
n = 1500

# cts X1 missing



source("SimData/GenData_normal_hm.R")

delta <- 2
rho <- 0

alpha = c(0,-0.5)
theta = rep(1,3)

if (delta == 0) {
  beta = c(-0.2, 0.5, 0.5, 0.5)
} else {
  beta = c(0.3, 0.5, 0.5, 0.5)
}

p.miss <- numeric(R)

seed = 1000
for ( r in 1:R ) {
  # seed <- round(runif(1,1,1e8))
  d <- GenData_normal_hm(N = n,
                         seed = seed + r,
                         alpha = alpha,
                         theta = theta,
                         beta = beta,
                         delta = delta,
                         rho = rho)
  p.miss[r] <- sum(is.na(d$Xp_obs[, 1]))/n
  file <- paste0("SimData/cts_rho", rho, "_delta", delta, sprintf("/data%04d",r))
  save(d,file=file)
}
mean(p.miss)

# bin X1 missing

source("SimData/GenData_binary_hm.R")

# delta <- 0.5
# rho <- 0

alpha = c(0,-0.5)
theta = rep(1,3)

if (delta == 0) {
  beta = c(-0.4, 0.5, 0.5, 0.5)
} else {
  beta = c(-1.2, 0.5, 0.5, 0.5)
}

# seed = 1000
for ( r in 1:R ) {
  # seed <- round(runif(1,1,1e8))
  d <- GenData_binary_hm(N = n,
                         seed = seed + r,
                         alpha = alpha,
                         theta = theta,
                         beta = beta,
                         delta = delta,
                         rho = rho)
  p.miss[r] <- sum(is.na(d$Xp_obs[, 1]))/n
  file <- paste0("SimData/bin_rho", rho, "_delta", delta, sprintf("/data%04d",r))
  save(d,file=file)
}
mean(p.miss)

# # cts X1, X2 missing
# 
# source("SimData/GenData_GeneralMissing_hm.R")
# 
# delta <- 0
# rho <- 0
# 
# alpha = c(0.2,-0.5)
# theta = rep(1,3)
# 
# if (delta == 0) {
#   beta = c(-0.2, 0.5, 0.5, 0.5)
#   beta2 <- c(0.1, 0.5, 0.5, 0.5)
# } else {
#   beta = c(-0.25, 0.5, 0.5, 0.5)
#   beta2 <- c(0.1, 0.5, 0.5, 0.5)
# }
# 
# 
# 
# seed = 1000
# for ( r in 1:R ) {
#   seed <- round(runif(1,1,1e8))
#   d <- GenData_GeneralMissing_hm(N = n,
#                                  seed = seed + r,
#                                  alpha = alpha,
#                                  theta = theta,
#                                  beta = beta,
#                                  beta2 = beta2,
#                                  delta = delta,
#                                  rho = rho)
#   sum(is.na(d$Xp_obs[, 1]))/n
#   sum(is.na(d$Xp_obs[, 2]))/n
#   file <- paste0("SimData/cts2_rho", rho, "_delta", delta, sprintf("/data%04d",r))
#   save(d,file=file)
# }

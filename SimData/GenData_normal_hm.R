## inputs
# N: sample size
# seed: reproducibility
# alpha: Used to generate Xp[, 1], the variable with missingness
## alpha to be deleted
# theta: y -- X association
# beta: missing -- X association
# sd.o: outcome's noise
# sd.x: missing variable value's noise

## intermediates
# p: dimension without intercept, including the variable with missingness
# Xp: N x p predictor matrix. The first variable is a linear combination of
#     other variables and has missingness.
# 

# 1000 data points were divided into 100 times 10 observations
# Then at most 9 out of 10 can be missing.


### MNAR
# X1: missing
# rX1: missing flag. 1 -- observed; 2 -- missing.
# alpha: X1 ~ 1 + X2
# theta: y ~ 1 + X1 + X2
# beta, delta: rX1 ~ beta(1, Y, X2, X3) + delta(X1)
# sd.o: outcome model noise
# sd.s: selection model noise
# sd.x: X2, X3 ~ N(0, sd.x^2)
# sd.y: Y's noise
# rho: covariance of sd.s and sd.o for selection
# _o: outcome
# _s: selection

GenData_normal_hm <- function(N,seed,alpha,theta,beta,delta,
                              sd.s=1,sd.o=1,sd.x=0.5,sd.y=1,
                              rho=0)
{
  set.seed(seed)
  
  p.y <- length(theta)
  p.s <- length(beta)
  p.o <- length(alpha)
  
  x2 <- rnorm(N, mean = 0, sd = sd.x)
  x3 <- rnorm(N, mean = 0, sd = sd.x)
  
  if (delta == 0) {
    covmat <- diag(2)
    covmat[1, 1] <- sd.s^2
    covmat[2, 2] <- sd.o^2
    covmat[1, 2] <- covmat[2, 1] <- rho * sd.o * sd.s
    noises <- mvtnorm::rmvnorm(N, sigma = covmat)
    
    x1 <- c(cbind(1, x2) %*% alpha + noises[, 2])
    
    muy <- c(cbind(1, x1, x2) %*% theta)
    y = muy + rnorm(N, sd = sd.y)
    
    miss.flag <- (cbind(1, y, x2, x3) %*% beta + noises[, 1] < 0)
  } else {
    x1 <- c(cbind(1, x2) %*% alpha) + rnorm(N, mean = 0, sd = sd.o)
    
    muy <- c(cbind(1, x1, x2) %*% theta)
    y = muy + rnorm(N, sd = sd.y)
    
    mu.obs <- c(cbind(1, y, x2, x3) %*% beta) + x1 * delta
    prob.missing <- 1/(1 + exp(mu.obs))
    miss.flag <- rbinom(N, 1, prob.missing)
  }
  
  #  pr = 1/(1+exp(-cbind(1,Xp) %*% theta))
  #  y = rbinom(N,1,pr)
  
  #missing data mechanism
  # mu.missing = cbind(1,y,x2,x3) %*% beta
  # prob.missing = matrix(1/(1+exp(-drop(mu.missing))),10)
  # miss.flag = matrix(rbinom(N,1,prob.missing),10)
  # sum.miss = apply(miss.flag,2,sum)
  # for ( i in 1:(N/10) )
  #   while ( sum.miss[i] > 9 )
  #   {
  #     miss.flag[,i] = rbinom(10,1,prob.missing[,i])
  #     sum.miss[i] = sum(miss.flag[,i])
  #   }
  
  Xp <- cbind(x1, x2, x3)
  Xp_obs = Xp
  Xp_obs[miss.flag==1,1] = NA
  # summary(Xp_obs[,1])
  
  d = new.env()
  d$model = "gaussian"
  d$N = N
  d$p.y = p.y
  d$p.s <- p.s
  d$p.o <- p.o
  d$y = y
  d$muy <- muy
  d$Xp = Xp
  d$Xp_obs = Xp_obs
  # d$prob.missing = prob.missing
  d$theta = theta
  d$beta = beta
  d$alpha = alpha
  d$rho <- rho
  d$delta <- delta
  
  # there is an additional item d$all added in lower level files SimXX.R
  
  d
}




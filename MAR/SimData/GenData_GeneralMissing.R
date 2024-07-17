
GenData_GeneralMissing <- function(N,seed,alpha,theta,beta1,beta2,beta3,sd.e=1,sd.x=1)
{
  set.seed(seed)
  
  p = length(theta) - 1
  
  Sigma = diag(1,3)
  Sigma[1,2] = Sigma[2,1] = Sigma[2,3] = Sigma[3,2] = Sigma[1,3] = Sigma[3,1] = 0.5
  cSigma = chol(Sigma)
  err = matrix(rnorm(3*N),N,3) %*% cSigma

  Xp_1 = matrix(rnorm((p-3)*N,sd=sd.x),N,p-3)
  Xp = cbind(drop(cbind(1,Xp_1)%*%alpha)+err,Xp_1)
  muy = cbind(1,Xp) %*% theta
  y = rnorm(N,muy,sd.e)
#  pr = 1/(1+exp(-cbind(1,Xp) %*% theta))
#  y = rbinom(N,1,pr)
  
  
  # missing data mechanism for X_1
  Xmiss = cbind(1,y,Xp_1)
  
  # missing for X_1
  mu.missing = Xmiss %*% beta1
  prob.missing = matrix(1/(1+exp(-drop(mu.missing))),10)
  miss.flag.1 = matrix(rbinom(N,1,prob.missing),10)
  sum.miss = apply(miss.flag.1,2,sum)
  for ( i in 1:(N/10) )
    while ( sum.miss[i] > 9 )
    {
      miss.flag.1[,i] = rbinom(5,1,prob.missing[,i])
      sum.miss[i] = sum(miss.flag.1[,i])
    }
  
  #missing for X_2
  mu.missing = Xmiss %*% beta2
  prob.missing = matrix(1/(1+exp(-drop(mu.missing))),10)
  miss.flag.2 = matrix(rbinom(N,1,prob.missing),10)
  sum.miss = apply(miss.flag.2,2,sum)
  for ( i in 1:(N/10) )
    while ( sum.miss[i] > 9 )
    {
      miss.flag.2[,i] = rbinom(5,1,prob.missing[,i])
      sum.miss[i] = sum(miss.flag.2[,i])
    }
  
  #missing for X_3
  mu.missing = Xmiss %*% beta3
  prob.missing = matrix(1/(1+exp(-drop(mu.missing))),10)
  miss.flag.3 = matrix(rbinom(N,1,prob.missing),10)
  sum.miss = apply(miss.flag.3,2,sum)
  for ( i in 1:(N/10) )
    while ( sum.miss[i] > 9 )
    {
      miss.flag.3[,i] = rbinom(5,1,prob.missing[,i])
      sum.miss[i] = sum(miss.flag.3[,i])
    }
  
  Xp_obs = Xp
  Xp_obs[miss.flag.1==1,1] = NA
  Xp_obs[miss.flag.2==1,2] = NA
  Xp_obs[miss.flag.3==1,3] = NA
  
  d = new.env()
  d$model = "gaussian"
  d$N = N
  d$p = p
  d$y = y
  d$Xp = Xp
  d$Xp_obs = Xp_obs
  d$alpha = alpha
  d$theta = theta
  d$beta = cbind(beta1,beta2,beta3)

  d
}



GenData_normal <- function(N,seed,alpha,theta,beta,sd.e=1,sd.x=1)
{
  set.seed(seed)
  
  p = length(theta) - 1

  Xp = matrix(runif(N*p,-3,3),N,p)
  Xp[,1] = rnorm(N,cbind(1,Xp[,-1])%*%alpha,sd.x)
  muy = cbind(1,Xp) %*% theta
  y = rnorm(N,muy,sd.e)
#  pr = 1/(1+exp(-cbind(1,Xp) %*% theta))
#  y = rbinom(N,1,pr)
  
  #missing data mechanism
  mu.missing = cbind(1,y,Xp[,-1]) %*% beta
  prob.missing = matrix(1/(1+exp(-drop(mu.missing))),10)
  miss.flag = matrix(rbinom(N,1,prob.missing),10)
  sum.miss = apply(miss.flag,2,sum)
  for ( i in 1:(N/10) )
    while ( sum.miss[i] > 9 )
    {
      miss.flag[,i] = rbinom(10,1,prob.missing[,i])
      sum.miss[i] = sum(miss.flag[,i])
    }

  Xp_obs = Xp
  Xp_obs[miss.flag==1,1] = NA
  
  d = new.env()
  d$model = "gaussian"
  d$N = N
  d$p = p
  d$y = y
  d$Xp = Xp
  d$Xp_obs = Xp_obs
  d$prob.missing = prob.missing
  d$theta = theta
  d$beta = beta
  
  d
}




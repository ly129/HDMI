
GenData_binary <- function(N,seed,alpha,theta,beta,sd.e=2)
{
  set.seed(seed)

  p = length(theta) - 1

  Xp_1 = matrix(runif(N*(p-1),-3,3),N,p-1)
  Xp = cbind(rbinom(N,1,prob=1/(1+exp(-cbind(1,Xp_1)%*%alpha))), Xp_1)
  muy = cbind(1,Xp) %*% theta
  y = rnorm(N,muy,sd.e)

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
  d$alpha = alpha
  d$theta = theta
  d$beta = beta
  
  d
}



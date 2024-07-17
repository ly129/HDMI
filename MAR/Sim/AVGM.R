source("Sim/PP.R")


AVGMLS = function(d,yidx,Xidx,n)
{
  p = length(Xidx)
  K = length(n)
  ni = diffinv(n)
  N = sum(n)
  
  fit = LS(d,yidx,Xidx,n)
  comm = 1
  
  df = sum(fit$df)
  SSE = sum(fit$SSE)
  beta = fit$beta%*%fit$df/df
  Cov = matrix(0,p,p)
  for ( k in 1:K )
    Cov = Cov + chol2inv(fit$cgram[,,k])*fit$df[k]^2/df^2
  cCov = chol(Cov)
  gram = chol2inv(cCov)
  cgram = chol(gram)

  list(beta=beta,df=df,SSE=SSE,gram=gram,cgram=cgram,comm=comm)
}


AVGMLogit = function(d,yidx,Xidx,n,beta0)
{
  p = length(Xidx)
  K = length(n)
  ni = diffinv(n)
  
  fit = Logit(d,yidx,Xidx,n,beta0, lam = 0.000) # lam = 0.005 for real data
  comm = 1
  
  df = sum(fit$df)
  beta = fit$beta%*%fit$df/df
  Cov = matrix(0,p,p)
  for ( k in 1:K )
    Cov = Cov + chol2inv(fit$cfisher[,,k])*fit$df[k]^2/df^2
  cCov = chol(Cov)
  fisher = chol2inv(cCov)
  cfisher = chol(fisher)

  list(beta=beta,fisher=fisher,cfisher=cfisher,comm=comm)
}

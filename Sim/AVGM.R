source("Sim/PP.R")

# d: data, of environment type
# yidx: input is midx (column index of X with missingness.
#       e.g. 3:5 or 3 depending on scenario.)
# Xidx: other variables in X, used to predict missingness in yidx.
## Xidx_s: column id of X_s, variables used to predict missing status
## Xidx_o: column id of X_o, variables used to predict missing variable value
# n: vector of length K, sample size of each site's data.
AVGMLS = function(d,yidx,Xidx_s,Xidx_o,n)
{
  p_s <- length(Xidx_s)
  p_o <- length(Xidx_o)
  p <- p_s + p_o + 2 # p (selection eq) + p-1 (outcome eq) + 2 (sigma & rho)
  K = length(n)
  ni = diffinv(n)
  N = sum(n)
  
  fit <- HMML_LS(d,yidx,Xidx_s,Xidx_o,n)
  # fit = LS(d,yidx,Xidx,n)
  comm = 1
  
  df = sum(fit$df)
  # SSE = sum(fit$SSE)
  beta = c(fit$beta%*%fit$df/df)
  Vb <- matrix(0,p,p)
  for (k in 1:K) {
    Vb <- Vb + fit$Vb[,,k] * fit$df[k]^2/df^2
  }
  # Cov = matrix(0,p,p)
  # for ( k in 1:K )
  #   Cov = Cov + chol2inv(fit$cgram[,,k])*fit$df[k]^2/df^2
  # cCov = chol(Cov)
  # gram = chol2inv(cCov)
  # cgram = chol(gram)

  list(beta=beta,Vb=Vb,df=df,comm=comm)
}

AVGMLogit = function(d,yidx,Xidx_s,Xidx_o,n)
{
  p_s <- length(Xidx_s)
  p_o <- length(Xidx_o)
  p <- p_s + p_o + 1 # p (selection eq) + p-1 (outcome eq) + 1 (rho) # no sigma
  K = length(n)
  ni = diffinv(n)
  N = sum(n)
  
  fit <- HMML_Logit(d,yidx,Xidx_s,Xidx_o,n)
  # fit = LS(d,yidx,Xidx,n)
  comm = 1
  
  df = sum(fit$df)
  # SSE = sum(fit$SSE)
  beta = c(fit$beta%*%fit$df/df)
  Vb <- matrix(0,p,p)
  for (k in 1:K) {
    Vb <- Vb + fit$Vb[,,k] * fit$df[k]^2/df^2
  }
  # Cov = matrix(0,p,p)
  # for ( k in 1:K )
  #   Cov = Cov + chol2inv(fit$cgram[,,k])*fit$df[k]^2/df^2
  # cCov = chol(Cov)
  # gram = chol2inv(cCov)
  # cgram = chol(gram)
  
  list(beta=beta,Vb=Vb,df=df,comm=comm)
}

# AVGMLogit = function(d,yidx,Xidx,n,beta0)
# {
#   p = length(Xidx)
#   K = length(n)
#   ni = diffinv(n)
#   
#   fit = Logit(d,yidx,Xidx,n,beta0)
#   comm = 1
#   
#   df = sum(fit$df)
#   beta = fit$beta%*%fit$df/df
#   Cov = matrix(0,p,p)
#   for ( k in 1:K )
#     Cov = Cov + chol2inv(fit$cfisher[,,k])*fit$df[k]^2/df^2
#   cCov = chol(Cov)
#   fisher = chol2inv(cCov)
#   cfisher = chol(fisher)
# 
#   list(beta=beta,fisher=fisher,cfisher=cfisher,comm=comm)
# }

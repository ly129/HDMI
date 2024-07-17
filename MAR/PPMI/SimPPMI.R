source("Sim/PP.R")

SimPPMI = function(scen,R,n,M,batch=0)
{
  K = length(n)
  method = sprintf("SIMI%d",K)
  
  if ( scen == 2 )
    MImethod = "logreg"
  else if ( scen == 3 )
    MImethod = rep("norm",3)
  else
    MImethod = "norm"
  
  if ( scen == 3 )
    midx = 3:5
  else
    midx = 3
  
  if ( scen == 3 )
    Theta = matrix(0,6,R)
  else
    Theta = matrix(0,3,R)
  
  comm = rep(0,R)
  
  for ( r in 1:R )
  {
    message(paste("Dataset=",r+batch,sep=""))
    load(sprintf("SimData/Scen%d/data%04d",scen,r+batch))
    
    fit = PPMI(d,M,midx,MImethod,n)
    Theta[,r] = fit$theta
    comm[r] = fit$comm
  }
  
  list(method=method,R=batch+1:R,scen=scen,Theta=Theta,comm=comm)
}


PPMI = function(d,M,midx,method,n)
{
  p = d$p+1
  q = length(midx)
  N = sum(n)
  K = length(n)
  ni = diffinv(n)
  thetas = matrix(0,p,M)
  
  d$all = cbind(d$y,1,d$Xp_obs)
  miss = is.na(d$all)
  
  comm = 0
  
  fit.imp = NULL
  for ( j in 1:q )
  {
    if ( method[j] == "norm" )
      fit.imp[[j]] = PPLS(d,midx[j],setdiff(1:(p+1),midx),n)
    else
      fit.imp[[j]] = PPLogit(d,midx[j],setdiff(1:(p+1),midx),n,c(0,d$alpha))
    comm = comm + fit.imp[[j]]$comm + 1
  }
  
  for ( m in 1:M )
  {
    for ( k in 1:K )
    {
      for ( j in 1:q )
      {
        if ( method[j] == "norm" )
        {
          sig = sqrt(1/rgamma(1,(fit.imp[[j]]$df+1)/2,(fit.imp[[j]]$SSE+1)/2))
          alpha = fit.imp[[j]]$beta + sig * backsolve(fit.imp[[j]]$cgram,rnorm(p))
          idx = ni[k] + which(miss[(ni[k]+1):ni[k+1],midx[j]])
          d$all[idx,midx[j]] = d$all[idx,setdiff(1:(p+1),midx)] %*% alpha + rnorm(length(idx),0,sig)
        }
        else
        {
          alpha = fit.imp[[j]]$beta + backsolve(fit.imp[[j]]$cfisher,rnorm(p))
          idx = ni[k] + which(miss[(ni[k]+1):ni[k+1],midx[j]])
          pr = 1 / (1 + exp(-d$all[idx,setdiff(1:(p+1),midx)] %*% alpha))
          d$all[idx,midx[j]] = rbinom(length(idx),1,pr)
        }          
      }
    }

    if ( d$model == "logistic" )
      fit = PPLogit(d,1,2:(p+1),n,d$theta)
    else
      fit = PPLS(d,1,2:(p+1),n)
    thetas[,m] = fit$beta
  }
  
  theta = apply(thetas,1,mean)

  list(theta=theta,comm=comm)
}

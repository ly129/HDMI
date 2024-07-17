source("Sim/PP.R")

SimIMI = function(scen,R,n,M,batch=0)
{
  K = length(n)
  method = sprintf("IMI%d",K)
  
  if ( scen == 2 ) {
    MImethod = "logreg"
  } else if ( scen == 3 ) {
    MImethod = rep("norm",3)
  } else {
    MImethod = "norm"
  }
  
  if ( scen == 3 ) {
    midx = 3:5
  } else {
    midx = 3
  }
  
  if ( scen == 3 ) {
    Theta = matrix(0,6,R)
  } else {
    Theta = matrix(0,3,R)
  }
  
  comm = rep(0,R)
  
  for ( r in 1:R ) {
    message(paste("Dataset=",r+batch,sep=""))
    if (scen == 1) {
      file <- paste0("../SimData/cts_rho", rho, "_delta",
                     delta, sprintf("/data%04d",r+batch))
    } 
    if (scen == 2) {
      file <- paste0("../SimData/bin_rho", rho, "_delta",
                     delta, sprintf("/data%04d",r+batch))
    }
    load(file)
    
    fit = IMI(d,M,midx,MImethod,n)
    Theta[,r] = fit$theta
    comm[r] = fit$comm
  }
  
  list(method=method,R=batch+1:R,scen=scen,Theta=Theta,comm=comm)
}



IMI = function(d,M,midx,method,n)
{
  # intercept, x1, x2 in analysis model, dimension p
  # intercept, y, x2, x3 in imputation model, dimension p + 1
  p = d$p.y
  q = length(midx)
  N = sum(n)
  K = length(n)
  ni = diffinv(n)
  thetas = matrix(0,p,M)
  
  d$all = cbind(d$y,1,d$Xp_obs)
  colnames(d$all)[1:2] <- c("y", "Intercept")
  miss = is.na(d$all)
  
  comm = 0
  
  fit.imp = vector(mode = "list", length = q)
  for ( j in 1:q )
  {
    if ( method[j] == "norm" )
      fit.imp[[j]] = LS(d,midx[j],setdiff(1:(p+2),midx),n)
    else
      fit.imp[[j]] = Logit(d,midx[j],setdiff(1:(p),midx),n,c(0,d$alpha), lam = 0.005)  # 1:p for real data, 1:p+2 for simulation, lam = 0.000 for simulation, lam = 0.005 for real data
  }
  
  for ( m in 1:M )
  {
    for ( k in 1:K )
    {
      for ( j in 1:q )
      {
        if ( method[j] == "norm" )
        {
          sig = sqrt(1/rgamma(1,(fit.imp[[j]]$df[k]+1)/2,(fit.imp[[j]]$SSE[k]+1)/2))
          alpha = fit.imp[[j]]$beta[,k] + sig * backsolve(fit.imp[[j]]$cgram[,,k],rnorm(p+1))
          idx = ni[k] + which(miss[(ni[k]+1):ni[k+1],midx[j]])
          d$all[idx,midx[j]] = d$all[idx,setdiff(1:(p+2),midx)] %*% alpha + rnorm(length(idx),0,sig)
        }
        else
        {
          alpha = fit.imp[[j]]$beta[,k] + backsolve(fit.imp[[j]]$cfisher[,,k],rnorm(p+1))
          idx = ni[k] + which(miss[(ni[k]+1):ni[k+1],midx[j]])
          pr = 1 / (1 + exp(-d$all[idx,setdiff(1:(p),midx)] %*% alpha)) # 1:p for real data, 1:p+2 for simulation
          d$all[idx,midx[j]] = rbinom(length(idx),1,pr)
        }          
      }
    }

    if ( d$model == "logistic" )
      fit = PPLogit(d,1,2:(p+1),sum(n),d$theta)
    else
      fit = PPLS(d,1,2:(p+1),sum(n))
    
    thetas[,m] = fit$beta
  }
  
  theta = apply(thetas,1,mean)
  
  list(theta=theta,comm=comm)
}



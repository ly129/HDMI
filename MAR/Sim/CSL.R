
CSLLS <- function(d,yidx,Xidx,n,ctr=1,lam=0.000)
{
  p = length(Xidx)
  K = length(n)
  N = sum(n)
  ni = diffinv(n)
  
  cc = apply(is.na(d$all[1:N,c(yidx,Xidx)]),1,sum) == 0
  d$cc = cc

  idx1 = ni[ctr] + which(cc[(ni[ctr]+1):ni[ctr+1]])
  ncc1 = length(idx1)
  X1 = matrix(d$all[idx1,Xidx],ncc1)
  y1 = d$all[idx1,yidx]
  XX1 = t(X1)%*%X1 + diag(lam,p)
  Xy1 = t(X1)%*%y1
  yy1 = sum(y1^2)
  cXX1 = chol(XX1)
  iXX1 = chol2inv(cXX1)
  beta1 = iXX1%*%Xy1

  idx = which(cc)
  ncc = length(idx)
  X = d$all[idx,Xidx]
  e = d$all[idx,yidx] - X%*%beta1
  Xe = t(X)%*%e

  beta = beta1 + iXX1%*%Xe*ncc1/ncc
  
  e1 = drop(y1 - X1%*%beta)
  SSE = sum(e1^2)*ncc/ncc1
  gram = XX1*ncc/ncc1
  cgram = chol(gram)

  if ( K == 1 )
    comm = 0
  else
    comm = 2
  
  list(beta=beta,df=ncc,SSE=SSE,gram=gram,cgram=cgram,comm=comm)
}



CSLLogit <- function(d,yidx,Xidx,n,beta0,ctr=1,lam=0.000,maxiter=50) # center set to which.max(n) = 1 for real data, lam = 0.000 for simulation, 0.005 for real data
{
  p = length(Xidx)
  K = length(n)
  N = sum(n)
  ni = diffinv(n)
  
  beta = beta0
  
  cc = apply(is.na(d$all[1:N,c(yidx,Xidx)]),1,sum) == 0
  d$cc = cc

  idx1 = ni[ctr] + which(cc[(ni[ctr]+1):ni[ctr+1]])
  ncc1 = length(idx1)
  X1 = d$all[idx1,Xidx]
  y1 = d$all[idx1,yidx]

  iter = 0
  while ( iter < maxiter )
  {
    iter = iter + 1
    
    xb = drop(X1%*%beta)
    pr = 1/(1+exp(-xb))
    H = t(X1)%*%(X1*pr*(1-pr)) + diag(lam,p)
    g = t(X1)%*%(y1-pr) - lam*beta
    Q = sum(y1*xb) + sum(log(1-pr[pr<0.5])) + sum(log(pr[pr>=0.5])-xb[pr>=0.5]) - lam*sum(beta^2)/2
    dir = chol2inv(chol(H))%*%g
    m = sum(dir*g)
    
    step = 1
    while (TRUE)
    {
      nbeta = beta + step*dir
      if ( max(abs(nbeta-beta)) < 1e-5 )
        break
      xb = drop(X1%*%nbeta)
      pr = 1/(1+exp(-xb))
      nQ = sum(y1*xb) + sum(log(1-pr[pr<0.5])) + sum(log(pr[pr>=0.5])-xb[pr>=0.5]) - lam*sum(nbeta^2)/2
      if ( nQ-Q > m*step/2 )
        break
      step = step / 2
    }
    
    if ( max(abs(nbeta-beta)) < 1e-5 )
      break
    beta = nbeta
  }
  
  idx = which(cc)
  ncc = length(idx)
  X = d$all[idx,Xidx]
  y = d$all[idx,yidx]
  pr = 1/(1+exp(-drop(X%*%beta)))
  g_off = - t(X)%*%(y-pr)*ncc1/ncc

  iter = 0
  while ( iter < maxiter )
  {
    iter = iter + 1
    
    xb = drop(X1%*%beta)
    pr = 1/(1+exp(-xb))
    H = t(X1)%*%(X1*pr*(1-pr)) + diag(lam,p)
    g = t(X1)%*%(y1-pr) - g_off - lam*beta
    Q = sum(y1*xb) + sum(log(1-pr[pr<0.5])) + sum(log(pr[pr>=0.5])-xb[pr>=0.5]) - sum(g_off*beta) - lam*sum(beta^2)/2
    dir = chol2inv(chol(H))%*%g
    m = sum(dir*g)
    
    step = 1
    while (TRUE)
    {
      nbeta = beta + step*dir
      if ( max(abs(nbeta-beta)) < 1e-5 )
        break
      xb = drop(X1%*%nbeta)
      pr = 1/(1+exp(-xb))
      nQ = sum(y1*xb) + sum(log(1-pr[pr<0.5])) + sum(log(pr[pr>=0.5])-xb[pr>=0.5]) - sum(g_off*nbeta) - lam*sum(nbeta^2)/2
      if ( nQ-Q > m*step/2 )
        break
      step = step / 2
    }
    
    if ( max(abs(nbeta-beta)) < 1e-5 )
      break
    beta = nbeta
  }
  beta = nbeta
  
  fisher = H * (ncc/ncc1)
  cfisher = chol(fisher)

  if ( K == 1 )
    comm = 0
  else
    comm = 2
  
  list(beta=beta,fisher=fisher,cfisher=cfisher,comm=comm)
}

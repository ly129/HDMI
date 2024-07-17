
LS <- function(d,yidx,Xidx,n,lam=0.000)
{
  p = length(Xidx)
  K = length(n)
  ni = diffinv(n)
  
  XX = array(0,c(p,p,K))
  df = rep(0,K)
  
  cA = array(0,c(p,p,K))
  beta = matrix(0,p,K)
  SSE = rep(0,K)
  
  cc = apply(is.na(d$all[,c(yidx,Xidx)]),1,sum) == 0
  d$cc = cc
  
  for ( k in 1:K )
  {
    idx = ni[k] + which(cc[(ni[k]+1):ni[k+1]])
    df[k] = length(idx)
    Xk = matrix(d$all[idx,Xidx],df[k],p)
    yk = matrix(d$all[idx,yidx],df[k],1)
    XX[,,k] = t(Xk)%*%Xk + diag(lam,p)
    Xy = t(Xk)%*%yk
    yy = sum(yk^2)

    cA[,,k] = chol(XX[,,k])
    beta[,k] = backsolve(cA[,,k],forwardsolve(t(cA[,,k]),Xy))
    SSE[k] = yy - sum(Xy*beta[,k])
  }

  comm = 0

  list(beta=beta,SSE=SSE,df=df,gram=XX,cgram=cA,comm=comm)
}


PPLS <- function(d,yidx,Xidx,n,lam=0.000)
{
  p = length(Xidx)
  K = length(n)
  N = sum(n)

  cc = apply(is.na(d$all[1:N,c(yidx,Xidx)]),1,sum) == 0
  d$cc = cc
  df = sum(cc)
  idx = which(cc)
  
  Xk = d$all[idx,Xidx]
  yk = d$all[idx,yidx]
  XX = t(Xk)%*%Xk + diag(lam,p)
  Xy = t(Xk)%*%yk
  yy = sum(yk^2)

  cA = chol(XX)
  beta = backsolve(cA,forwardsolve(t(cA),Xy))
  SSE = yy - sum(Xy*beta)

  if ( K == 1 )
    comm = 0
  else
    comm = 1
  
  list(beta=beta,SSE=SSE,df=df,gram=XX,cgram=cA,comm=comm)
}


Logit = function(d,yidx,Xidx,n,beta0,lam=0.000,maxiter=50)
{
  p = length(Xidx)
  K = length(n)
  ni = diffinv(n)
  
  H = array(0,c(p,p,K))
  cH = array(0,c(p,p,K))
  beta = matrix(beta0,p,K)
  
  cc = apply(is.na(d$all[,c(yidx,Xidx)]),1,sum) == 0
  d$cc = cc
  df = rep(0,K)
  
  for ( k in 1:K )
  {
    idx = ni[k] + which(cc[(ni[k]+1):ni[k+1]])
    ncc = length(idx)
    df[k] = ncc
    Xk = matrix(d$all[idx,Xidx],ncc,p)
    yk = matrix(d$all[idx,yidx],ncc,1)
    
    iter = 0
    while ( iter < maxiter )
    {
      iter = iter + 1
      
      xb = drop(Xk%*%beta[,k])
      prk = 1/(1+exp(-xb))
      wk = prk*(1-prk)
      H[,,k] = t(Xk)%*%(Xk*wk) + diag(lam,p)
      g = t(Xk)%*%(yk-prk) - lam*beta[,k]
      Q = sum(yk*xb) + sum(log(1-prk[prk<0.5])) + sum(log(prk[prk>=0.5])-xb[prk>=0.5]) - lam*sum(beta[,k]^2)/2
      cH[,,k] = chol(H[,,k])
      dir = backsolve(cH[,,k],forwardsolve(t(cH[,,k]),g))
      m = sum(dir*g)
      
      step = 1
      while (TRUE)
      {
        nbeta = beta[,k] + step*dir
        if ( max(abs(nbeta-beta[,k])) < 1e-5 )
          break
        xb = drop(Xk%*%nbeta)
        prk = 1/(1+exp(-xb))
        nQ = sum(yk*xb) + sum(log(1-prk[prk<0.5])) + sum(log(prk[prk>=0.5])-xb[prk>=0.5]) - lam*sum(nbeta^2)/2
        if ( nQ-Q > m*step/2 )
          break
        step = step / 2
      }
      
      if ( max(abs(nbeta-beta[,k])) < 1e-5 )
        break
      beta[,k] = nbeta
    }
    beta[,k] = nbeta
  }
  
  comm = 0
  
  list(beta=beta,df=df,fisher=H,cfisher=cH,comm=comm)
}


PPLogit = function(d,yidx,Xidx,n,beta0,lam=0.000,maxiter=50)
{
  p = length(Xidx)
  K = length(n)
  N = sum(n)
  
  beta = beta0
  
  cc = apply(is.na(d$all[1:N,c(yidx,Xidx)]),1,sum) == 0
  d$cc = cc
  idx = which(cc)
  ncc = length(idx)
  
  X = d$all[idx,Xidx]
  y = d$all[idx,yidx]
  
  comm = 0
  while ( comm < maxiter )
  {
    comm = comm + 1
    
    Xb = drop(X%*%beta)
    pr = 1/(1+exp(-Xb))
    w = pr*(1-pr)
    H = t(X)%*%(X*w) + diag(lam,p)
    g = t(X)%*%(y-pr) - lam*beta
    Q = sum(y*Xb) + sum(log(1-pr[pr<0.5])) + sum((log(pr[pr>=0.5])-Xb[pr>=0.5])) - lam*sum(beta^2)/2
    
    cH = chol(H)
    dir = backsolve(cH,forwardsolve(t(cH),g))
    m = sum(dir*g)
    
    step = 1
    while (TRUE)
    {
      nbeta = beta + step*dir
      if ( max(abs(nbeta-beta)) < 1e-5 )
        break
      Xb = drop(X%*%nbeta)
      pr = 1/(1+exp(-Xb))
      nQ = sum(y*Xb) + sum(log(1-pr[pr<0.5])) + sum((log(pr[pr>=0.5])-Xb[pr>=0.5])) - lam*sum(nbeta^2)/2
      if ( nQ-Q > m*step/2 )
        break
      step = step / 2
    }
    
    if ( max(abs(nbeta-beta)) < 1e-5 )
      break
    beta = nbeta
  }
  beta = nbeta
  
  if ( K == 1 )
    comm = 0
  
  list(beta=beta,df=ncc,fisher=H,cfisher=cH,comm=2*comm)
}

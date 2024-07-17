## intermediate
# XX: X transpose X from K sites.
# df: 
# lam: added to the diagonal of XtX, for positive semi-definiteness
# yidx -- column id of outcome y in d$all
# Xidx -- column ids of predictors X in d$all
# note: in d$all, 1st column is y, 2nd is intercept, rest is X, see SimXX.R

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

# ### HM2_LS: Heckman's two-step estimation for continuous missing variable
# ### input
# # d: environment that includes the data
# # yidx: column id of y (i.e. the variable with missing)
# # Xidx_s: column id of X_s, variables used to predict missingness (selection)
# # Xidx_o: column id of X_o, variables used to predict outcome
# ## Galimard 2018 says that X_s and X_o cannot be exactly the same
# ## as this may produce collinearity issues and possibly erroneous estimation
# ## to avoid this, include at least one supplementary variable in the selection
# # n: sample sizes of the data sources
# HM2_LS <- function(d,yidx,Xidx_s,Xidx_o,n)
# {
#   p = length(Xidx)
#   K = length(n)
#   ni = diffinv(n)
#   
#   XX = array(0,c(p,p,K))
#   df = rep(0,K)
#   
#   cA = array(0,c(p,p,K))
#   # original beta has dimension p, the addition dim is IMR in Heckman 2-stage
#   beta = matrix(0,p + 1,K)
#   SSE = rep(0,K)
#   
#   cc = apply(is.na(d$all[,c(yidx,Xidx)]),1,sum) == 0
#   d$cc = cc
#   
#   for ( k in 1:K )
#   {
#     idk = (ni[k] + 1):ni[k+1]
#     dallk <- as.data.frame(d$all[idk, ])
#     var_sel <- paste0("V", Xidx_s)
#     var_out <- paste0("V", Xidx_o)
#     selec <- paste0("V", yidx)
#     dallk$selec <- !is.na(dallk[, selec])
#     
#     formula_s <- as.formula(paste("selec ~",
#                                   paste(var_sel, collapse= "+"), " - 1"))
#     formula_o <- as.formula(paste(selec, "~",
#                                   paste(var_out, collapse= "+"), " - 1"))
#     
#     hm2fit <- sampleSelection::heckit2fit(selection = formula_s,
#                                           outcome = formula_o,
#                                           data = dallk)
#     
#     
#   }
#   
#   comm = 0
#   
#   list(beta=beta,SSE=SSE,df=df,gram=XX,cgram=cA,comm=comm)
# }

### HMML_LS: Heckman's (one-step) maximum likelihood estimation
### for continuous missing variable
# Inputs:
## Xidx_s: column id of X_s, variables used to predict missing status
## Xidx_o: column id of X_o, variables used to predict missing variable value
# output:
## beta: fitted coefficients for imputation (theta.star is the estimate of rho)
## Vb: variance-covariance matrix of fitted beta
HMML_LS <- function(d,yidx,Xidx_s,Xidx_o,n)
{
  K = length(n)
  ni = diffinv(n)
  p_s <- length(Xidx_s)
  p_o <- length(Xidx_o)
  p <- p_s + p_o + 2 # p (selection eq) + p-1 (outcome eq) + 2 (sigma & rho)
  
  beta = matrix(0,p,K)
  Vb <- array(dim = c(p, p, K))
  df <- integer(K)
  
  for ( k in 1:K )
  {
    idk = (ni[k] + 1):ni[k+1]
    dallk <- as.data.frame(d$all[idk, ])
    
    var_sel <- names(dallk)[Xidx_s]
    var_out <- names(dallk)[Xidx_o]
    dallk$select <- !is.na(dallk[, yidx])
    
    formula_s <- as.formula(paste("select ~",
                                  paste(var_sel, collapse= "+"), " - 1"))
    formula_o <- as.formula(paste(names(dallk)[yidx], "~",
                                  paste(var_out, collapse= "+"), " - 1"))
    
    hmmlfit <- GJRM::copulaSampleSel(list(formula_s, formula_o),
                                     data = dallk, fp = TRUE)
    beta[, k] <- hmmlfit$coefficients
    Vb[, , k] <- hmmlfit$Vb
    # gradb[, k] <- hmmlfit$fit$gradient # for CSL
    df[k] <- hmmlfit$n.sel
  }
  
  comm = 0
  
  return(list(beta=beta,df=df,Vb=Vb,comm=comm))
}

### HMML_Logit: Heckman's (one-step) maximum likelihood estimation
### for binary missing variable
# Inputs are the same as in HMML_LS
HMML_Logit <- function(d,yidx,Xidx_s,Xidx_o,n)
{
  K = length(n)
  ni = diffinv(n)
  p_s <- length(Xidx_s)
  p_o <- length(Xidx_o)
  p <- p_s + p_o + 1 # p (selection eq) + p-1 (outcome eq) + 1 (rho) # no sigma for binary
  
  beta = matrix(0,p,K)
  Vb <- array(dim = c(p, p, K))
  df <- integer(K)
  
  for ( k in 1:K )
  {
    idk = (ni[k] + 1):ni[k+1]
    dallk <- as.data.frame(d$all[idk, ])
    
    var_sel <- names(dallk)[Xidx_s]
    var_out <- names(dallk)[Xidx_o]
    dallk$select <- !is.na(dallk[, yidx])
    
    formula_s <- as.formula(paste("select ~",
                                  paste(var_sel, collapse= "+"), " - 1"))
    formula_o <- as.formula(paste(names(dallk)[yidx], "~",
                                  paste(var_out, collapse= "+"), " - 1"))
    
    hmmlfit <- GJRM::SemiParBIV(list(formula_s, formula_o),
                                data = dallk, Model = "BSS",
                                fp = TRUE)
    beta[, k] <- hmmlfit$coefficients
    Vb[, , k] <- hmmlfit$Vb
    df[k] <- hmmlfit$n.sel
  }
  
  comm = 0
  
  return(list(beta=beta,df=df,Vb=Vb,comm=comm))
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
  
  # # debug, check XX singularity
  # cat("XX = ", XX, "\n")
  
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

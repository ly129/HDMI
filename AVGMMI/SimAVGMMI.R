source("Sim/AVGM.R")

SimAVGMMI = function(scen,R,n,M,batch=0)
{
  K = length(n)
  method = sprintf("AVGMMI%d",K)
  
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
    # yidx <- 3
    Xidx_s <- c(1, 2, 4, 5)
    Xidx_o <- c(1, 2, 4)
  }
  
  if ( scen == 3 ) {
    Theta = matrix(0,6,R)
  } else {
    Theta = matrix(0,3,R)
  }
  comm = rep(0,R)
  
  for ( r in 1:R )
  {
    message(paste("Dataset=",r+batch,sep=""))
    if (scen == 1) {
      file <- paste0("SimData/cts_rho", rho, "_delta", delta, sprintf("/data%04d",r+batch))
    } 
    if (scen == 2) {
      file <- paste0("SimData/bin_rho", rho, "_delta", delta, sprintf("/data%04d",r+batch))
    }
    load(file)
    
    fit = AVGMMI(d,M,midx,Xidx_s,Xidx_o,MImethod,n)
    Theta[,r] = fit$theta
    comm[r] = fit$comm
  }
  
  list(method=method,R=batch+1:R,scen=scen,Theta=Theta,comm=comm)
}


# midx: column index of X with missingness. e.g. 3:5 or 3 depending on scenario.
# M = 20 imputed datasets
# j = 1...q, number of missing variables
# k = 1... K sites
AVGMMI = function(d,M,midx,Xidx_s,Xidx_o,method,n)
{
  p_s <- length(Xidx_s)
  p_o <- length(Xidx_o)
  q = length(midx)
  N = sum(n)
  K = length(n)
  ni = diffinv(n)
  thetas = matrix(0,d$p.y,M)
  
  d$all = cbind(d$y,1,d$Xp_obs)
  colnames(d$all)[1:2] <- c("y", "Intercept")
  miss = is.na(d$all)
  
  comm = 0
  
  fit.imp = vector(mode = "list", length = q)
  for ( j in 1:q ) {
    if ( method[j] == "norm" )
      fit.imp[[j]] = AVGMLS(d,midx[j],Xidx_s,Xidx_o,n)
    else
      fit.imp[[j]] = AVGMLogit(d,midx[j],Xidx_s,Xidx_o,n)
    comm = comm + fit.imp[[j]]$comm + 1
  }
  
  for ( m in 1:M ) {
    for ( k in 1:K ) {
      idx = ni[k] + which(miss[(ni[k]+1):ni[k+1],midx[j]])
      if (length(idx) != 0) {
        for ( j in 1:q ) {
          if ( method[j] == "norm" ) {
            alpha <- mvtnorm::rmvnorm(n = 1,
                                      mean = fit.imp[[j]]$beta,
                                      sigma = fit.imp[[j]]$Vb)
            idx = ni[k] + which(miss[(ni[k]+1):ni[k+1],midx[j]])
            var_sel <- d$all[idx, Xidx_s]
            betaxstar_sel <- c(var_sel %*% alpha[1:p_s])
            
            var_out <- d$all[idx, Xidx_o]
            betaxstar_out <- c(var_out %*% alpha[(1+p_s):(p_s + p_o)])
            
            rho <- alpha[p_s + p_o + 2]
            rho <- pmax(pmin(rho, 100), -100)
            rho_star <- (exp(2 * rho) - 1)/(1 + exp(2 * rho))
            
            sigma <- alpha[p_s + p_o + 1]
            sigma_star <- exp(sigma)
            
            y.star <- betaxstar_out +
              sigma_star*rho_star*
              (-dnorm(betaxstar_sel)/(pnorm(-betaxstar_sel))) +
              rnorm(length(idx),0, sigma_star)
            d$all[idx,midx[j]] <- y.star
          } else {
            alpha <- mvtnorm::rmvnorm(n = 1,
                                      mean = fit.imp[[j]]$beta,
                                      sigma = fit.imp[[j]]$Vb)
            idx = ni[k] + which(miss[(ni[k]+1):ni[k+1],midx[j]])
            var_sel <- d$all[idx, Xidx_s]
            betaxstar_sel <- c(var_sel %*% alpha[1:p_s])
            
            var_out <- d$all[idx, Xidx_o]
            betaxstar_out <- c(var_out %*% alpha[(1+p_s):(p_s + p_o)])
            
            rho <- alpha[p_s + p_o + 1]
            rho <- pmax(pmin(rho, 100), -100)
            rho_star <- (exp(2 * rho) - 1)/(1 + exp(2 * rho))
            
            p.star <- pbivnorm::pbivnorm(betaxstar_out, -betaxstar_sel, -rho_star)/pnorm(-betaxstar_sel)
            p.star <- pmax(pmin(p.star, 1-1e-6), 1e-6)
            
            y.star <- rbinom(length(idx), 1, p.star)
            d$all[idx,midx[j]] <- y.star
          }          
        }
      }
    }
    
    if ( d$model == "logistic" ) {
      fit = PPLogit(d,1,2:(p_o+1),n,d$theta)
    } else {
      fit = PPLS(d,1,2:(p_o+1),n)
    }
    
    thetas[,m] = fit$beta
  }
  
  theta = apply(thetas,1,mean)
  
  list(theta=theta,comm=comm)
}



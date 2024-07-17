source("Sim/PP.R")

SimGS <- function(scen,R,n,batch=0)
{
  method = "CD"
  
  if ( scen == 3 ) {
    Theta = matrix(0,6,R)
  } else {
    Theta = matrix(0,3,R)
  }
  comm = rep(0,R)
  
  for ( r in 1:R ) {
    message(paste("Dataset=",r+batch,sep=""))
    if (scen == 1) {
      file <- paste0("SimData/cts_rho", rho, "_delta", delta, sprintf("/data%04d",r+batch))
    } 
    if (scen == 2) {
      file <- paste0("SimData/bin_rho", rho, "_delta", delta, sprintf("/data%04d",r+batch))
    }
    if (scen == 3) {
      file <- paste0("SimData/cts2_rho", rho, "_delta", delta, sprintf("/data%04d",r+batch))
    }
    load(file)
    
    ### y, intercept, x, a total of p+2 columns
    d$all = cbind(d$y,1,d$Xp)
    colnames(d$all)[1:2] <- c("y", "Intercept")
    
    if ( d$model == "logistic" )
      fit = PPLogit(d,1,2:(d$p.y+1),n,d$theta)
    else
      fit = PPLS(d,1,2:(d$p.y+1),n)
    Theta[,r] = fit$beta
    comm[r] = fit$comm
  }
  
  list(method=method,R=batch+1:R,scen=scen,Theta=Theta,comm=comm)
}

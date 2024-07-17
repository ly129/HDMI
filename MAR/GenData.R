


R = 1000 # number of MC data
n = 1000

##################
# Scenario 1 & 2 #
##################

source("SimData/GenData_normal.R")
alpha = c(0.2,-0.5)
theta = rep(1,3)
beta = c(0.3,-0.2,0.1) 
seed = 1000
for ( r in 1:R )
{
  d = GenData_normal(n,seed+r,alpha,theta,beta)
  save(d,file=sprintf("SimData/Scen1/data%04d",r))
}



##################
# Scenario 3 & 4 #
##################

source("SimData/GenData_binary.R")
alpha = c(0.2,-0.5)
theta = rep(1,3)
beta = c(0.3,-0.2,0.1) 
seed = 2000
for ( r in 1:R )
{
  d = GenData_binary(n,seed+r,alpha,theta,beta)
  save(d,file=sprintf("SimData/Scen2/data%04d",r))
}



##################
# Scenario 5 & 6 #
##################

source("SimData/GenData_GeneralMissing.R")
alpha = c(0.3,-0.3,-0.1)
theta = rep(1,6)
beta1 = c(-1.0,-0.4,-0.1,-0.2)   
beta2 = c(-0.8,-0.6,0.2,0.4)  
beta3 = c(-0.8,-1,0.4,0.3)   
seed = 3000
for ( r in 1:R )
{
  d = GenData_GeneralMissing(n,seed+r,alpha,theta,beta1,beta2,beta3)
  save(d,file=sprintf("SimData/Scen3/data%04d",r))
}





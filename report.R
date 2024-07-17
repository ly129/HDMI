##########
# Report #
##########


report <- function(res)
{
  L = length(res)
  # cat(" & & & method & bias & sd & rmse & \\#c\\\\\n")
  cat("Method\tBias\tSD\tRMSE\tComm\n")
  
  for ( l in 1:L )
  {
    p = dim(res[[l]]$Theta)[1]
    theta = rep(1,p)
    # theta <- c(0,1,1)
    R = length(res[[l]]$R)
    m = apply(res[[l]]$Theta,1,mean)
    bias = sqrt(sum((m - theta)^2))
    se = sqrt(sum((res[[l]]$Theta-m)^2)/R)
    rmse = sqrt(sum((res[[l]]$Theta-theta)^2)/R)
    comm = mean(res[[l]]$comm)
    #cat(sprintf(" & & & %s & & %.3f & %.3f & %.3f & & %.1f\\\\\n",
    cat(sprintf("%s\t%.3f\t%.3f\t%.3f\t%.1f\n",
                res[[l]]$method,bias,se,rmse,comm))
  }
}
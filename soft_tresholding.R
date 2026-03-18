soft_tresholding <- function(Y, c)
{
  y = abs(Y) - c
  ind = which(0 < y)
  z = rep(0, length(Y)) 
  z[ind] <- y[ind]
  z = sign(Y)*z
  return(z)
}


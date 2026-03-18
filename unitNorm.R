unitNorm = function(M){
  M <- as.matrix(M)
  if (ncol(M)>1){
    ssm <- apply(M, 2, function(x) sqrt(sum(x^2)))
    Ms <- M %*% diag(1/ssm)
  }
  else{
    ssm <- sqrt(sum(M^2))
    Ms <- M * (1/ssm)
  }
  return(list(Ms = Ms, ssm = ssm))
}

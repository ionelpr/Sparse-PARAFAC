Phi = function(a, b){
  require(MASS)
  if ((is.matrix(a)) & (is.matrix(b))) {
    p = ginv(diag(diag(t(a) %*% a), nrow = ncol(a))^0.5) %*% 
      t(a) %*% b %*% ginv(diag(diag(t(b) %*% b), nrow = ncol(b))^0.5)
  }
  else {
    p = ginv((t(a) %*% a)^0.5) %*% t(a) %*% b %*% ginv((t(b) %*% 
                                                            b)^0.5)
  }
  return(p)
}

phi = function(a, b){
  require(MASS)
  if ((is.matrix(a)) & (is.matrix(b))) {
    p = ginv(diag(diag(t(a) %*% a), nrow = ncol(a))^0.5) %*% 
      t(a) %*% b %*% ginv(diag(diag(t(b) %*% b), nrow = ncol(b))^0.5)
  }
  else {
    p = ginv((t(a) %*% a)^0.5) %*% t(a) %*% b %*% ginv((t(b) %*% 
                                                          b)^0.5)
  }
  return(p)
}


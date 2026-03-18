randMat = function(r, c){
  mat <- matrix((rbinom(r*c, 1, .5) * 2-1) * runif(r*c, 0.5, 1), r, c)
  return(mat)
}


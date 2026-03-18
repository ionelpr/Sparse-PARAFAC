require(ThreeWay)
orthNorm <- function(n,r){
  return(orth(matrix(runif(n * r, 0, 1), nrow = n) - 0.5))
}
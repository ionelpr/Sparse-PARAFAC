ed <- function(s) {
  ev <- eigen(s)
  values <- ev$values
  vectors <- ev$vectors
  
  #idx <- order(values, decreasing = TRUE)
  #values <- values[idx]
  #vectors <- vectors[, idx]
  
  vectors <- nrm(vectors)
  
  list(u = vectors, v = diag(values))
}
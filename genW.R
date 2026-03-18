genW = function(n, k, s){ # n = rows; k = columns; s = simplicity, degree of overlap between components.
  # partition
  cluster <- sample(c(1:k, sample(1:k, n-k, replace = T)))
  W <- diag(k)[cluster,]
  
  # overlap
  overlap <- round(s * (n*k - n), 0)
  zeros <- which(W == 0, arr.ind = TRUE)
  nz <- nrow(zeros)
  zeros <- zeros[sample(1:nz), ]
  
  #fill
  pos <- sample(1:nz, overlap)
  W[zeros[pos, ]] <- 1
  
  return(W) 
}


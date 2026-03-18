# Parafac with L0 penalty on A, B, C

source("ssq.R")
require(L0Learn)
require(ThreeWay)

CPL0LAM <- function(X, n, m, p, r, lamA, lamB, lamC, conv, maxit, start, A = NULL, B = NULL, C = NULL){
  #  min || X - A (C kr B)' || + lamA ||A||_0 + lamB ||B||_0 + lamC ||C||_0 
  #  (kr: Khatri-Rao product)
  # Input: 
  #        X: data (supermatrix n x mp)
  #        n,m,p: order
  #        r: dimensionality of solution
  #        lamA, lamB, lamC: penalty weights for A, B, C
  #   	   conv: convergence criterion, e.g., 1e-6
  #        maxit: maximal number of iterations
  # Output:
  #        A,B,C: component matrices
  #        f: loss function value
  #		     fp: fit percentage
  #	       cputime: computation time
  # Require:
  #        ThreeWay
  #        L0Learn
  cputime <- system.time({
    ssx <- ssq(X)
    # start by SVD
    if(start == 0){
      A <- eigen(X%*%t(X))$vectors[, 1:r] 
      Z <- permnew(X, n, m, p)
      B <- eigen(Z%*%t(Z))$vectors[, 1:r]
      Z <- permnew(Z, m, p, n)
      C <- eigen(Z%*%t(Z))$vectors[, 1:r] 
    }
    # random start
    if(start == 1){
      A <- orth(matrix(runif(n*r), n, r) - 0.5)
      B <- orth(matrix(runif(m*r), m, r) - 0.5) 
      C <- orth(matrix(runif(p*r), p, r) - 0.5) 
    }
    
    H <- matrix(0, r, r*r) 
    for(ii in 1:r){
      H[ii, (ii-1)*r + ii] <- 1
    }
    H1 <- permnew(H, r, r, r)
    H1 <- permnew(B%*%H1, m, r, r)
    H1 <- permnew(C%*%H1, p, r, m)
    
    f <- ssq(X - A%*%H1) + lamA*sum(A!=0) + lamB*sum(B!=0) + lamC*sum(C!=0)
    fold <- f + 2*conv*f
    iter <- 0
    BB <- t(B)%*%B
    CC <- t(C)%*%C
    
    while((fold - f > conv*f | iter < 2) & f > conv^2 & iter < maxit){
      fold <- f
      # update A
      y_A <- as.matrix(c(X))
      X_A <- t(H1) %x% diag(n)
      updA <- coef(L0Learn.fit(X_A, y_A, lambdaGrid = list(lamA), maxSuppSize = n*r, atol = 1e-3, rtol = 1e-3, intercept = FALSE))
      A <- matrix(updA, n, r)
      if (any(colSums(A==0) == n))
      {
        print("At least one column of A has all 0's: reduce lamA")
        A <- matrix(NA, n, r)
        B <- matrix(NA, m, r)
        C <- matrix(NA, p, r)
        f <- NA
        ftiter <- NA
        break
      }
      AA <- t(A)%*%A
      # Update B
      Z <- permnew(X, n, m, p)
      y_B <- as.matrix(c(Z))
      H1 <- permnew(H, r, r, r)
      H1 <- permnew(C %*% H1, p, r, r)
      H1 <- permnew(A %*% H1, n, r, p)
      X_B <- t(H1) %x% diag(m)
      updB <- coef(L0Learn.fit(X_B, y_B, lambdaGrid = list(lamB), maxSuppSize = m*r, atol = 1e-3, rtol = 1e-3, intercept = FALSE))
      B <- matrix(updB, m, r)
      if (any(colSums(B==0) == m))
      {
        print("At least one column of B has all 0's: reduce lamB")
        A <- matrix(NA, n, r)
        B <- matrix(NA, m, r)
        C <- matrix(NA, p, r)
        f <- NA
        ftiter <- NA
        break
      }
      BB <- t(B) %*% B
      # Update C
      Z <- permnew( permnew(X, n, m, p), m, p, n )
      y_C <- as.matrix(c(Z))
      H1 <- permnew(H, r, r, r)
      H1 <- permnew(A %*% H1, n, r, r)
      H1 <- permnew(B %*% H1, m, r, n)
      X_C <- t(H1) %x% diag(p)
      updC <- coef(L0Learn.fit(X_C, y_C, lambdaGrid = list(lamC), maxSuppSize = p*r, atol = 1e-3, rtol = 1e-3, intercept = FALSE))
      C <- matrix(updC, p, r)
      if (any(colSums(C==0) == p))
      {
        print("At least one column of C has all 0's: reduce lamC")
        A <- matrix(NA, n, r)
        B <- matrix(NA, m, r)
        C <- matrix(NA, p, r)
        f <- NA
        ftiter <- NA
        break
      }
      CC <- t(C) %*% C
      # Evaluate
      H1 <- permnew(H, r, r, r)
      H1 <- permnew(B %*% H1, m, r, r)
      H1 <- permnew(C %*% H1, p, r, m)
      f <- ssq(X - A %*% H1) + lamA*sum(A!=0) + lamB*sum(B!=0) + lamC*sum(C!=0)
      iter <- iter + 1
    }
  })
  tripcos = min(Phi(A, A) * Phi(B, B) * Phi(C, C))
  fp <- 100 - 100*ssq(X - A%*%H1)/(ssq(X))
  return(list(A = A, B = B, C = C, f = f, fp = fp, cputime = cputime, lamA = lamA, lamB = lamB, lamC = lamC, tripcos = tripcos))  
}

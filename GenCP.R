source("randMat.R")
source("unitNorm.R")
source("genW.R")

SIMPmat <- rbind(
  c(1-0.25, 1-0.25, 1-0.25), # (same small)
  c(1-0.50, 1-0.50, 1-0.50), # (same large)
  c(1-0.50, 1-0.25, 1-0.25), # (different - large A)
  c(1-0.25, 1-0.50, 1-0.25), # ( different - large B)
  c(1-0.25, 1-0.25, 1-0.50), # ( different - large C)
  c(1-0.50, 1-0.50, 1-0.25), # ( different - large A & B) 
  c(1-0.50, 1-0.25, 1-0.50), # ( different  - large A & C)
  c(1-0.25, 1-0.50, 1-0.50)) # ( different  - large B & C) 

IJKmat <- rbind(c(15,15,15), c(30, 15, 15), c(30, 30, 15),
                c(30, 30, 30), c(60, 15, 15), c(60, 30, 15), c(60, 60, 15),
                c(60, 30, 30), c(60, 60, 30), c(60, 60, 60))

GenCP = function(grow){
  
  # gridrow structure: "M", "IJK", "S", "SIMP", "snr", "seed"
  
  I <- IJKmat[grow$IJK, 1]
  J <- IJKmat[grow$IJK, 2]
  K <- IJKmat[grow$IJK, 3]
  simp <- SIMPmat[grow$SIMP,]
  set.seed(grow$seed)
  
  # plain loadings
  A <- randMat(I, grow$S)
  B <- randMat(J, grow$S)
  C <- randMat(K, grow$S)
  
  # Generate weight matrices (rows, column, simplicity a.k.a. components overlap)
  Wa <- genW(I, grow$S, simp[1])
  Wb <- genW(J, grow$S, simp[2])
  Wc <- genW(K, grow$S, simp[3])
  
  # Sparsify the matrices according to W matrices
  
  AW = A*Wa
  BW = B*Wb
  CW = C*Wc
  
  Va <- diag(colSums(AW^2)); Va[row(Va) != col(Va)] <- 0 # 
  Vb <- diag(colSums(BW^2)); Vb[row(Vb) != col(Vb)] <- 0 # 
  Vc <- diag(colSums(CW^2)); Vc[row(Vc) != col(Vc)] <- 0 # 
  
  # identity core
  G = matrix(0, grow$S, grow$S*grow$S)
  for(ii in 1:grow$S){
    G[ii, (ii-1)*grow$S + ii] <- (1/sqrt(Va[ii,ii])) * (1/sqrt(Vb[ii,ii])) * (1/sqrt(Vc[ii,ii]))
  }
  
  # reconstruct the pure signal
  X <- AW %*% G %*% t(kronecker(CW, BW))
  
  # build the error (~ ensuring the snr)
  E <- matrix(rnorm(I*J*K), I, J*K)
  epsilon <- sqrt(sum(X^2) / (sum(E^2) * grow$snr))
  
  # add the noise to the signal
  Xa <- X + epsilon * E
  
  # count the zeros on each loading matrix
  pA <- sum(AW==0)
  pB <- sum(BW==0)
  pC <- sum(CW==0)
  
  return(list(Xa = Xa, Atrue = AW, Btrue = BW, Ctrue = CW, I = I, J = J, K = K, S = grow$S, pA = pA, pB = pB, pC = pC, Wa = Wa, Wb = Wb, Wc = Wc, Va = Va, Vb = Vb, Vc = Vc))
}


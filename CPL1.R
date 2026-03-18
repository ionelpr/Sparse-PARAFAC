source("PTD_L1L1L1.R")
source("L1_update_unconstrained.R")
source("inits.R")
source("soft_tresholding.R")
source("product.R")
source("unitNorm.R")
require(rTensor)
require(ThreeWay)

# X:            dato in forma tensoriale (oggetto rTensor); 
# a/b/c_init:   inizializzazioni di A, B, C (comp x oggetti); 
# c1/c2/c3:     intensita delle penalty; 
# Niter:        nr. iterazioni

CPL1 = function(X, n, m, p, r, lamA, lamB, lamC, maxit = 500, conv = 1e-6, start=1, A = NULL, B = NULL, C=NULL) # funzione principale
{
cputime = system.time({
  
  if(start == 1){
    a_init <- t(unitNorm(matrix(inits(n*r), n, r))$Ms)
    b_init <- t(unitNorm(matrix(inits(m*r), m, r))$Ms)
    c_init <- t(unitNorm(matrix(inits(p*r), p, r))$Ms)
  }

  if(start == 2){
    if(is.null(A) || is.null(B) || is.null(C)){
      stop("A, B, C must be provided when start = 2.")
      }
    else{
      a_init = t(A)
      b_init = t(B)
      c_init = t(C)
      }
  }
  
  if(start == 0){
    if (n >= r) {
      AUT = eigen(X %*% t(X))
      a_init = t(AUT$vectors[, 1:r])
    }
    else {
      A = orth(matrix(runif(r * r, 0, 1), nrow = r) - 
                 0.5)
      a_init = t(A[1:n, ])
    }
    Z = permnew(X, n, m, p)
    if (m >= r) {
      AUT = eigen(Z %*% t(Z))
      b_init = t(AUT$vectors[, 1:r])
    }
    else {
      B = orth(matrix(runif(r * r, 0, 1), nrow = r) - 
                 0.5)
      b_init = t(B[1:m, ])
    }
    Z = permnew(Z, m, p, n)
    if (p >= r) {
      AUT = eigen(Z %*% t(Z))
      c_init = t(AUT$vectors[, 1:r])
    }
    else {
      C = orth(matrix(runif(r * r, 0, 1), nrow = r) - 
                 0.5)
      c_init = (C[1:p, ])
    }
  }
  X <- k_fold(X, m = 1, modes = c(n, m, p))
  
  # dato A, B, C, ricavo w (lambda nell'articolo): w_k = <X, a[,k] o b[,k] o c[,k]>
  w_hat = rep(0, times = r)
  
  for(j in 1:r)
  {
    w_hat[j] = product(X, c(1,2,3), a_init[j,], b_init[j,], c_init[j,]) 
  }
  
  # dato A, B, C, w, ricavo Y (ricostruzione di X) Y = sum_k w_k * a[,k] o b[,k] o c[,k]
  Y = as.tensor(array(0,dim(X)))
  
  obj_value_prev =  fnorm(Y - X) + lamA*sum(abs(a_init)) + lamB*sum(abs(b_init)) + lamC*sum(abs(c_init))
  
  for(j in 1:r)
  {
    lizt <- list('mat' = as.matrix(a_init[j,]),'mat2' = as.matrix(b_init[j,]),'mat3'= as.matrix(c_init[j,]))
    aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
    
    Y = Y + w_hat[j]*aux 
  }
  
  for(iter in 1:maxit)
  {  
    for(j in 1:r)
    {
      # 
      lizt <- list('mat' = as.matrix(a_init[j,]),'mat2' = as.matrix(b_init[j,]),'mat3'= as.matrix(c_init[j,]))
      aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
      
      Y  = Y -  w_hat[j]*aux # Ricostruzione totale - ricostruzione via k-esima componente
      Z =  X - Y # dato residuo: dato - ricostruzione via componenti 1:k 
      
      temp2 = PTD_L1L1L1(Z = Z, c1 = lamA, c2 = lamB, c3 = lamC, Niter = 5) # ricavo a[,k], b[,k], c[,k] ottimali per il dato residuo
      
      # ottengo la componente k-esima w_k * a[,k] o b[,k] o c[,k] che meglio approssima il dato residuo Z
      a_init[j,] = as.vector(temp2$u)
      b_init[j,] = as.vector(temp2$v)      
      c_init[j,] = as.vector(temp2$w)
      w_hat[j] = product(Z, c(1,2,3), a_init[j,], b_init[j,], c_init[j,])
      
      lizt <- list('mat' = as.matrix(a_init[j,]), 'mat2' = as.matrix(b_init[j,]), 'mat3'= as.matrix(c_init[j,]))
      aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
      
      # aggiorno la k-esima componente della ricostruzione Y
      Y = Y + w_hat[j]*aux 
    }
    
    obj_value =  fnorm(Y - X)^2 + lamA*sum(abs(a_init)) + lamB*sum(abs(b_init)) + lamC*sum(abs(c_init))
    fp = ( fnorm(Y)^2 / fnorm(X)^2 ) * 100 
    tripcos = min(Phi(t(a_init), t(a_init)) * Phi(t(b_init), t(b_init)) * Phi(t(c_init), t(c_init)))
    if( (obj_value - obj_value_prev) / (obj_value_prev +.00001) < conv )
    { break }
    
    obj_value_prev = obj_value
    
  }
})
  return(list(A = t(a_init), B = t(b_init), C = t(c_init), w_hat = w_hat, lamA = lamA, lamB = lamB, lamC = lamC, loss = obj_value, f = obj_value, fp = fp, cputime = cputime[1], tripcos = tripcos))
}


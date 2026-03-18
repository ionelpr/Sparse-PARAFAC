require(rTensor)
source("L1_update_unconstrained.R")

PTD_L1L1L1 <- function(Z,c1,c2,c3,Niter) 
{
  # inizializzazione della componente a, b, c (per approssimare Z)
  u = inits(dim(Z)[1])
  v = inits(dim(Z)[2])
  w = inits(dim(Z)[3])
  
  u = u / sqrt(sum(u*u))   
  v = v / sqrt(sum(v*v))
  w = w / sqrt(sum(w*w))
  
  # warm-up
  
  for(j in 1:5)
  {
    lizt <- list('mat2' = v,'mat3' =w)
    u =   ttl(Z, lizt, ms = c(2,3)) 
    u = t(as.matrix(u@data))
    u =  u/norm(u,"F")
    
    lizt <- list('mat2' = u,'mat3' =w)
    v =   ttl(Z, lizt, ms = c(1,3)) 
    v = t(as.matrix(v@data))
    v =  v/norm(v,"F")
    
    lizt <- list('mat2' = u,'mat3' =v)
    w =   ttl(Z, lizt, ms = c(1,2)) 
    w = t(as.matrix(w@data))
    w =  w/norm(w,"F")
  }
  
  u_prev = u; v_prev = v; w_prev = w;
  
  lizt <- list('mat1' = u_prev, 'mat2' = v_prev, 'mat3' = w_prev)
  fprev =  ttl(Z, lizt, ms = c(1,2,3)) 
  f = fprev
  # la funzione obiettivo è da massimizzare: f/fprev misura la similitudine (?) tra la ricostruzione (lizt) e Z
  fprev = fprev@data 
  
  # iterativamente, aggiorno la componente in a, b, c con una regressione penalizzata
  for( iter in 1:Niter)
  {
    
    u = t(as.matrix(L1_update_unconstrained(Z, v, w, c1, 1))) 
    if(length(which(u!=0))==0)
    { iter = 42; break }
    
    v = t(as.matrix(L1_update_unconstrained(Z, u, w, c2, 2)))
    if(length(which(v!=0))==0)
    { iter = 42; break }
    
    w = t(as.matrix(L1_update_unconstrained(Z, u, v, c3, 3))) 
    if(length(which(w!=0))==0)
    { iter = 42; break }
    
    lizt <- list('mat1' = u, 'mat2' = v, 'mat3' = w)
    f =  ttl(Z, lizt, ms = c(1,2,3)) 
    f = f@data 
    
    aux = abs( (f-fprev) / (fprev+0.0001) )
    
    if(aux < 1e-6)
    { break }
    
    fprev = f
    u_prev = u; v_prev = v; w_prev = w;
  }
  
  u = t(as.matrix(as.vector(u)))
  v = t(as.matrix(as.vector(v)))
  w = t(as.matrix(as.vector(w)))
  
  lizt <- list('mat1' = u,'mat2' = v,'mat3' = w)
  d = ttl(Z, lizt, ms = c(1,2,3))  
  
  return(list(iter = iter, u = u, v = v, w = w, d = d, f = f)) 
}


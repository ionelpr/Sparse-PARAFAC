require(rTensor)

product = function(X,ind,u,v,w)
{
  u = t(as.matrix(u))
  v = t(as.matrix(v))
  w = t(as.matrix(w))
  
  if(length(ind ) ==2)
  {
    lizt <- list('mat2' = u,'mat3' =v)
    u =   ttl(X, lizt, ms = ind)
    return( as.vector(u@data))
  }
  
  lizt <- list('mat' = u,'mat2' =v,'mat3'=w)
  u =   ttl(X, lizt, ms = ind)
  
  return( as.vector(u@data))
}
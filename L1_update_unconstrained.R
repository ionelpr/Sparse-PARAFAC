# regressione penalizzata
require(rTensor)
source("soft_tresholding.R")

L1_update_unconstrained = function(X, u, v, c3, index)
{  
  if(index == 1)
  {   
    lizt <- list('mat1' = u,'mat2' = v)
    
    Y = ttl(X, lizt, ms = c(2,3))
    Y = unfold(Y, row_idx = 1, col_idx = c(2,3))
    Y = t(Y@data)
  }
  
  if(index == 2)
  {   
    lizt <- list('mat1' = u,'mat2' = v)
    
    Y = ttl(X, lizt, ms = c(1,3))
    Y = unfold(Y, row_idx = 2, col_idx = c(1,3))
    Y = t(Y@data)
  }
  
  if(index == 3)
  {   
    lizt <- list('mat1' = u,'mat2' = v)
    
    Y = ttl(X, lizt, ms = c(1,2))
    Y = unfold(Y, row_idx = 3, col_idx = c(1,2))
    Y = t(Y@data)
  }
  
  aux = soft_tresholding(Y, c3)
  
  if( norm(aux,"F")==0)
  {
    return(matrix(0, dim(aux)[1], dim(aux)[2]))
  }
  
  return(t(aux / norm(aux,"F")))
}  


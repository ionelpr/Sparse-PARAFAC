# Simulations: Comparison of CP-ALS / CP-L1 / CP-L0 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(doParallel)

source("gold.R")
source("Phi.R")

# Design Matrix: 4 Design Variables. M dataset for each cell of the design matrix.

## Variable 1. data size
IJKmat <- rbind(c(15,15,15), c(30, 15, 15), c(30, 30, 15),
                c(30, 30, 30), c(60, 15, 15), c(60, 30, 15), c(60, 60, 15),
                c(60, 30, 30), c(60, 60, 30), c(60, 60, 60))
IJK <- c(1:nrow(IJKmat))

## Variable 2. components
S <- c(2, 3, 4)

## Variable 3. degree of overlap: The higher s, the lower the number of zeros. 
#     The number of zeros per loading matrix is: n * (k-1) - (round(s*n*(k-1), 0))
SIMPmat <- rbind(
c(1-0.25, 1-0.25, 1-0.25), # (same small)
c(1-0.50, 1-0.50, 1-0.50), # (same large)
c(1-0.50, 1-0.25, 1-0.25), # (different - large A)
c(1-0.25, 1-0.50, 1-0.25), # ( different - large B)
c(1-0.25, 1-0.25, 1-0.50), # ( different - large C)
c(1-0.50, 1-0.50, 1-0.25), # ( different - large A & B) 
c(1-0.50, 1-0.25, 1-0.50), # ( different  - large A & C)
c(1-0.25, 1-0.50, 1-0.50)) # ( different  - large B & C) 
SIMP <- c(1:nrow(SIMPmat))

## Variable 4. noise as signal-to-noise ratio
snr <- c(0.25, 0.5, 1)

## M replicas
M <- c(1:10)

# Create a Grid
Grid <- expand.grid(M, IJK, S, SIMP, snr)
names(Grid) <- c("M", "IJK", "S", "SIMP", "snr")

seed <- 1:nrow(Grid)
Grid <- cbind(Grid,seed)

N <- nrow(Grid)
cores <- detectCores()-2
cl <- makeCluster(cores, outfile = "") # 
registerDoParallel(cl)

### incremental saving
step <- cores
s <- 1

models_out <- list()
totSteps <- ceiling(N/step) 

for(i in 1:totSteps){ #  ?
  n = min(N, s + step - 1)
  cat(s:n)
  
  models_out[s:n] <- foreach(j = s:n, .packages = c("ThreeWay", "MASS", "glmnet", "L0Learn", "rTensor")) %dopar% {
    gold(Grid[j, ])
  }
  
  save(models_out, file = "thesims_out/out.RData")
  
  s = s + step
}

stopCluster(cl)

# from list to matrix
models_mat <- do.call(rbind, lapply(models_out, unlist))
save(models_mat, file="thesims_out/out_mat.RData")





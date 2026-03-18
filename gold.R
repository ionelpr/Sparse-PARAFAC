source("genCP.R")
source("CPs.R")
source("Eval.R")

gold = function(grow){
  
  set.seed(grow$seed)
  
  # generate one cell of the design matrix 
  gencp <- GenCP(grow)
  
  # Inserire tryCatch
  
  # apply the CP models
  cps <- CPs(gencp = gencp, grow = grow)
  
  # Evaluate ALS, L1, L0, L0rat
  EvalALS = Eval(gencp = gencp, cpsModBest = cps$ALS$best)
  
  EvalL1 = Eval(gencp = gencp, cpsModBest = cps$L1$best)
  EvalL1rat = Eval(gencp = gencp, cpsModBest = cps$L1rat)
  
  EvalL0 = Eval(gencp = gencp, cpsModBest = cps$L0$best)
  EvalL0rat = Eval(gencp = gencp, cpsModBest = cps$L0rat)
  
  # vectorize outputs
  
  base <- c(grow, gencp$pA, gencp$pB, gencp$pC)

  ALS <- c(unlist(EvalALS), cps$ALS$best$fp, cps$ALS$locott, cps$ALS$best$cputime, cps$ALS$time, cps$ALS$best$tripcos)

  L1 <- c(unlist(EvalL1), cps$L1$best$fp, cps$L1$locott, cps$L1$best$cputime, cps$L1$time, cps$L1$best$lamA, cps$L1$best$lamB, cps$L1$best$lamC, cps$L1$best$tripcos, sum(cps$L1$best$A==0), sum(cps$L1$best$B==0), sum(cps$L1$best$C==0), cps$L1$best$NRDimSel, cps$L1$best$err_ds)
  L1rat <- c(unlist(EvalL1rat), cps$L1rat$fp, cps$L1rat$cputime[1], cps$L1rat$lambda, cps$L1rat$tripcos)
  
  L0 <- c(unlist(EvalL0), cps$L0$best$fp, cps$L0$locott, cps$L0$best$cputime[1], cps$L0$time, cps$L0$best$lamA, cps$L0$best$lamB, cps$L0$best$lamC, cps$L0$best$tripcos, sum(cps$L0$best$A==0), sum(cps$L0$best$B==0), sum(cps$L0$best$C==0), cps$L0$best$NRDimSel, cps$L0$best$err_ds)
  L0rat <- c(unlist(EvalL0rat), cps$L0rat$fp, cps$L0rat$cputime[1], cps$L0rat$lambda, cps$L0rat$tripcos)
  # naming
  
  names(ALS)[(10+3):(14+3)] <- c("fp", "PercLocOtt", "tBest", "tN", "tripcos")
  
  names(L1)[(10+3):(22+3)] <- c("fp", "PercLocOtt", "tBest", "tN", "lamA", "lamB", "lamC", "tripcos", "zA", "zB", "zC", "nrDSleq2", "err_ds")
  names(L1rat)[(10+3):(15+3)] <- c("fp", "tBest", "lamA", "lamB", "lamC", "tripcos")
  
  names(L0)[(10+3):(22+3)] <- c("fp", "PercLocOtt", "tBest", "tN", "lamA", "lamB", "lamC", "tripcos",  "zA", "zB", "zC", "nrDSleq2", "err_ds")
  names(L0rat)[(10+3):(15+3)] <- c("fp", "tBest", "lamA", "lamB", "lamC", "tripcos")
  
  
  names(base)[7:9] <- c("pA", "pB", "pC")
  
  names(ALS) <- paste0("ALS_", names(ALS))
  
  names(L1) <- paste0("L1_", names(L1))
  names(L1rat) <- paste0("L1rat_", names(L1rat))
  
  names(L0) <- paste0("L0_", names(L0))
  names(L0rat) <- paste0("L0rat_", names(L0rat))
  
  return(c(base, ALS, L1, L1rat, L0, L0rat))
}


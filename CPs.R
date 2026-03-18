source("DimSelectorPenalty.R")
source("CpLamPath.R")
source("CPL1.R")
source("CPL0LAM.R")
source("orthNorm.R")

require(ThreeWay)

CPs <- function(gencp, grow){
  # per maggiore leggibilità del codice
  I <- gencp$I; J <- gencp$J; K <- gencp$K; S <- grow$S
  Xa <- gencp$Xa; pA <- gencp$pA; pB <- gencp$pB; pC = gencp$pC
  nrep <- 20 # random starts
  nlambda <- 50

  ### >>> CP-ALS
  set.seed(grow$seed)
  #ALStime <- system.time({
    ALS <- lapply(1:nrep, function(x) {
      Astart <- orthNorm(I, S); Bstart <- orthNorm(J, S); Cstart <- orthNorm(K,S);
      CPfunc(X = Xa, n = I, m = J, p = K, ort1 = 1, ort2 = 1, ort3 = 1, r = S, start = 2, A = Astart, B = Bstart, C = Cstart, conv = 1e-6, maxit = 500)
    })
  #})
  
  fps <- unlist((lapply(ALS, function(x) x$fp)))
  
  ALStime <- mean(unlist((lapply(ALS, function(x) x$cputime))))
  
  locottALS <- (sum(fps <= (max(fps) * 0.9999)) / nrep) * 100
  
  ALSbest <- ALS[[which.max(fps)]]
  
  VabcALS <- rbind(colSums(ALSbest$A^2), colSums(ALSbest$B^2), colSums(ALSbest$C^2))
  ALSbest$Vabc <- VabcALS
  
  ALSout <- list(best = ALSbest, locott = locottALS, time = ALStime)
  
  ### >>> CP-L1
  # outt <- matrix(0, nrow = 0, ncol = 9)
  # colnames(outt) <- c("n0.A", "n0.B", "n0.C", "Fit", "S", "nfp", "lamA", "lamB", "lamC")
  # for(s in (S-1):(S+1)){
  #   ou <- cp_lam_path(X = Xa, n = I, m = J, p = K, r = s, L0 = F, K = nlambda, max_backtrack = 500, t_init = 0.1, t_decay = 0.75, backtrack_factor = 0.95, verbose = T)$CPRF
  #   outt <- rbind(outt, ou)
  #}
  clp <- cp_lam_path(X = Xa, n = I, m = J, p = K, r = S, L0 = F, K = nlambda, max_backtrack = 500, t_init = 1, t_decay = 0.9, backtrack_factor = 0.9, verbose = T)
  outt <- clp$CPRF
  
  flag = 0 
  err_ds <- 0
  set.seed(grow$seed)
  while(flag==0){
    out <- DimSelectorPenalty(outt); l1 = as.numeric(out$st_max[7]); l2 = as.numeric(out$st_max[8]); l3 = as.numeric(out$st_max[9]); 
  
      L1 <- Filter(Negate(is.null),lapply(1:nrep, function(x) {
        tryCatch({
          withCallingHandlers({
            Astart <- orthNorm(I, S); Bstart <- orthNorm(J, S); Cstart <- orthNorm(K,S);
            CPL1(X = Xa, n = I, m = J, p = K, S, lamA = l1, lamB = l2, lamC = l3, conv = 1e-6, maxit = 500, start = 2, A = Astart, B = Bstart, C = Cstart)
          }, warning = function(w) invokeRestart("muffleWarning"))
        },
        error = function(e) NULL)
      })
      )
    
    fps <- unlist((lapply(L1, function(x) x$fp)))
    if(is.null(fps)) {
      idx <- which(apply(outt, 1, function(x) all(x == out$st_max[-10])))
      outt <- outt[-idx,]
      err_ds <- 1
    }
    else{ flag = 1 }
  }
  
  L1time <- mean(unlist((lapply(L1, function(x) x$cputime))), na.rm = T)
  flag = nrep - length(fps)
  if(flag != 0){ fps <- c(fps, rep(-Inf, flag)) }
  locottL1 <- (sum(fps <= (max(fps) * 0.9999)) / nrep) * 100
  
  L1best <- L1[[which.max(fps)]]
  if(nrow(out$out_st)<3) { L1best$NRDimSel = 1} else {L1best$NRDimSel = 0}
  L1best$err_ds <- err_ds
  
  VabcL1 <- rbind(colSums(L1best$A^2), colSums(L1best$B^2), colSums(L1best$C^2))
  L1best$Vabc <- VabcL1
  
  L1out <- list(best = L1best, locott = locottL1, time = L1time)
  
  # >>> rational L1
  set.seed(grow$seed)
  #L1rat <- CPL1(X = Xa, n = I, m = J, p = K, r = S, lamA = l1, lamB = l2, lamC = l3, A = ALSbest$A, B = ALSbest$B, C = ALSbest$C, conv = 1e-6, maxit = 500)
  target <- c(l1, l2, l3)
  # Take, from the solutions within the path, the one that has the lambdas used for L1 in the previous case
  # use that as rational  
  idx <- which(sapply(clp$solutions, function(sol) {
    lam <- sol$lambda
    !is.null(lam) && length(lam) == length(target) && all(lam == target)
    }))
  L1rat <- clp$solutions[[idx]]
  VabcL1rat <- rbind(colSums(L1rat$A^2), colSums(L1rat$B^2), colSums(L1rat$C^2))
  L1rat$Vabc <- VabcL1rat
  
  ### >>> CP-L0
  # outt <- matrix(0, nrow = 0, ncol = 9)
  # colnames(outt) <- c("n0.A", "n0.B", "n0.C", "Fit", "S", "nfp", "lamA", "lamB", "lamC")
  # for(s in (S-1):(S+1)){
  #   ou <- cp_lam_path(X = Xa, n = I, m = J, p = K, r = s, L0 = T, K = nlambda, max_backtrack = 500, t_init = 0.1, t_decay = 0.75, backtrack_factor = 0.95, verbose = T)$CPRF
  #   outt <- rbind(outt, ou)
  #}
  
  clp <- cp_lam_path(X = Xa, n = I, m = J, p = K, r = S, L0 = T, K = nlambda, max_backtrack = 500, t_init = 1, t_decay = 0.9, backtrack_factor = 0.9, verbose = T, sA = J*K, sB = I*K, sC = I*J)
  outt <- clp$CPRF
  
  flag = 0
  err_ds <- 0
  set.seed(grow$seed)
  while(flag==0){
    out <- DimSelectorPenalty(outt); ##
    lamA = as.numeric(out$st_max[7]); lamB = as.numeric(out$st_max[8]); lamC = as.numeric(out$st_max[9]); 
    
    L0time <- system.time({
      L0 <- Filter(Negate(is.null),lapply(1:nrep, function(x) {
                     tryCatch({
                       withCallingHandlers({
                         Astart <- orthNorm(I, S); Bstart <- orthNorm(J, S); Cstart <- orthNorm(K,S);
                         CPL0LAM(Xa, I, J, K, S, lamA, lamB, lamC, conv = 1e-6, maxit = 500, start = 2, A = Astart, B = Bstart, C = Cstart)
                       }, warning = function(w) invokeRestart("muffleWarning"))
                     },
                     error = function(e) NULL)
                                  })
                  )
            })
    
    fps <- unlist((lapply(L0, function(x) x$fp))) 
    if(is.null(fps)) {
      idx <- which(apply(outt, 1, function(x) all(x == out$st_max[-10])))
      outt <- outt[-idx,]
      err_ds <- 1
    }
    else{ flag = 1 }
  }
  L0time <- mean(unlist((lapply(L0, function(x) x$cputime))), na.rm = T)
  flag = nrep - length(fps)
  if(flag != 0){ fps <- c(fps, rep(-Inf, flag)) }
  locottL0 <- (sum(fps <= (max(fps) * 0.9999)) / nrep) * 100
  
  L0best <- L0[[which.max(fps)]]
  if(nrow(out$out_st)<3) { L0best$NRDimSel = 1} else {L0best$NRDimSel = 0}
  L0best$err_ds <- err_ds
  VabcL0 <- rbind(colSums(L0best$A^2), colSums(L0best$B^2), colSums(L0best$C^2))
  L0best$Vabc <- VabcL0
  L0out <- list(best = L0best, locott = locottL0, time = L0time)
  
  # >>> rational L0
  #set.seed(grow$seed)
  #L0rat <- CPL0LAM(X = Xa, n = I, m = J, p = K, r = S, lamA = lamA, lamB = lamB, lamC = lamC, start = 2, A = ALSbest$A, B = ALSbest$B, C = ALSbest$C, conv = 1e-6, maxit = 500)
  target <- c(lamA, lamB, lamC)
  idx <- which(sapply(clp$solutions, function(sol) {
    lam <- sol$lambda
    !is.null(lam) && length(lam) == length(target) && all(lam == target)
  }))
  L0rat <- clp$solutions[[idx]]
  
  VabcL0rat <- rbind(colSums(L0rat$A^2), colSums(L0rat$B^2), colSums(L0rat$C^2))
  L0rat$Vabc <- VabcL0rat
  
  return(list(ALS = ALSout, L1 = L1out, L1rat = L1rat, L0 = L0out, L0rat = L0rat))
}


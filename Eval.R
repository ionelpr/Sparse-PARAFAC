require(ThreeWay)

Eval = function(gencp, cpsModBest){
  ### Error of loading-matrix reconstruction
  
  A <- gencp$Atrue; B <- gencp$Btrue; C <- gencp$Ctrue;
  I <- gencp$I; J <- gencp$J; K <- gencp$K; S <- gencp$S;
  Wa <- gencp$Wa; Wb <- gencp$Wb; Wc <- gencp$Wc;
  Aest <- cpsModBest$A; Best <- cpsModBest$B; Cest <- cpsModBest$C;
  Vabc <- cpsModBest$Vabc # i vettori diagonali della norma al quadrato delle componenti sono impilati riga per riga in Vabc: 3 x S
  
  Wa.est = (Aest != 0)*1
  Wb.est = (Best != 0)*1
  Wc.est = (Cest != 0)*1
  
  perm <- perms(S)
  nperm <- nrow(perm)
  rec <- rep(NA, nperm)
  
  for (np in 1:nperm){
    p <- perm[np,]
    rec[np] <- sum( abs(diag(t(Aest[,p]) %*% A)) / sqrt( diag(gencp$Va) * Vabc[1,p] ) ) / S
  } 
  
  pm <- perm[which.max(rec),]
  
  signs = c()
  for(s in 1:S){
    
    signs[s] = sum( diag( t(Aest[,pm[s]]) %*% A[,s]  / sqrt( gencp$Va[s,s] * Vabc[1,pm[s]] ) ))
    if(signs[s]<0) { Aest[,pm[s]] = -Aest[,pm[s]]; Best[,pm[s]] = -Best[,pm[s]] }
  }
  
  for(s in 1:S){
    
    signs[s] = sum( diag( t(Best[,pm[s]]) %*% B[,s]  / sqrt( gencp$Vb[s,s] * Vabc[2,pm[s]] ) ))
    if(signs[s]<0) {Best[,pm[s]] = -Best[,pm[s]]; Cest[,pm[s]] = -Cest[,pm[s]]}
    
  }
  # check su segno su ciascuna colonna di Aest[,pm] (ciclo for)
  # [ sum( (diag(t(Aest[,p]) %*% A)) / sqrt( diag(gencp$Va) * Vabc[1,p] ) ) / S   ><  sum( (diag(t((-1)*Aest[,p]) %*% A)) / sqrt( diag(gencp$Va) * Vabc[1,p] ) ) / S]
  # se colonna s di Aest[,pm] deve essere moltiplicata *-1 (cioè se disuglianza vale con < )
  # allora colonna s di Best[,pm] deve essere moltiplicata *-1
  
  # check su segno su ciascuna colonna di Best[,pm] (ciclo for)
  # se colonna s di Best[,pm] deve essere moltiplicata *-1
  # allora colonna s di Cest[,pm] deve essere moltiplicata *-1
  
  # Tucker's Congruence Coefficient (on un-normalized matrices)
  recA.max <- sum((diag(t(Aest[,pm]) %*% A)) / sqrt( diag(gencp$Va) * Vabc[1,pm] ) ) / S
  recB.max <- sum((diag(t(Best[,pm]) %*% B)) / sqrt( diag(gencp$Vb) * Vabc[2,pm] ) ) / S
  recC.max <- sum((diag(t(Cest[,pm]) %*% C)) / sqrt( diag(gencp$Vc) * Vabc[3,pm] ) ) / S

  # NUOVI OK
  # recA.max <- sum( (diag(t(Aest[,pm]) %*% A)) / sqrt( diag(gencp$Va) * Vabc[1,pm] ) ) / S
  # recB.max <- sum( (diag(t(Best[,pm]) %*% B)) / sqrt( diag(gencp$Vb) * Vabc[2,pm] ) ) / S
  # recC.max <- sum( (diag(t(Cest[,pm]) %*% C)) / sqrt( diag(gencp$Vc) * Vabc[3,pm] ) ) / S
  
  # Correctly Estimated Zeros
  # Using the permutation according to the best rec
  cezA <- (sum(Wa == 0 & Wa.est[,pm] == 0) / sum(Wa==0)) * 100
  cezB <- (sum(Wb == 0 & Wb.est[,pm] == 0) / sum(Wb==0)) * 100
  cezC <- (sum(Wc == 0 & Wc.est[,pm] == 0) / sum(Wc==0)) * 100
  
  # True non-zero estimated as zero.
  cnzA <- mean(abs(unitNorm(A)$M - unitNorm(Aest[,pm])$M)[Wa==1])
  cnzB <- mean(abs(unitNorm(B)$M - unitNorm(Best[,pm])$M)[Wb==1])
  cnzC <- mean(abs(unitNorm(C)$M - unitNorm(Cest[,pm])$M)[Wc==1])
  
  # Mean Absolute Deviation on the positions of (true) zeros (MADPZ)
  
  ########################## 
  # madpzA <- mean(abs(unitNorm(A)$M - unitNorm(Aest[,pm])$M)[Wa==0])
  #########################
  madpzA <- mean(abs(unitNorm(A)$M - unitNorm(Aest[,pm])$M)[Wa==0])
  madpzB <- mean(abs(unitNorm(B)$M - unitNorm(Best[,pm])$M)[Wb==0])
  madpzC <- mean(abs(unitNorm(C)$M - unitNorm(Cest[,pm])$M)[Wc==0])
  
  out <- list(cezA = cezA, cezB = cezB, cezC = cezC, tccA = recA.max, tccB = recB.max, tccC = recC.max, madpzA = madpzA, madpzB = madpzB, madpzC = madpzC, cnzA = cnzA, cnzB = cnzB, cnzC = cnzC)
  return(out)
}


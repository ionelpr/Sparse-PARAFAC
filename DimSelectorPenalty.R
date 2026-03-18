require(ThreeWay)

DimSelectorPenalty <- function(out) 
{
  
  tmp1 = nrow(out)
  fpi = min(out[, 6])		# soluzione piu' parsimoniosa (minor nfp)
  
  # in fpi metto tutti i possibili valori distinti di nfp (dal minore al maggiore)
  # for (i in 1:tmp1) {
  #   if (out[i, 6] != max(fpi)) {
  #     fpi = c(fpi, out[i, 6])
  #   }
  # }
  
  fpi = sort(unique(out[,6])) # 1
  
  # in out_best vado a mettere tutte le possibili soluzioni di interesse
  # cioe', dato nfp, tengo soltanto la soluzione a cui corrisponde il massimo Fit
  
   tmp3 = length(fpi)
   out_best = matrix(0, tmp3, ncol(out))
   for (i in 1:tmp3) {
	 for (j in 1:tmp1) {
	 	if (fpi[i] == out[j, 6]) {
	 		if (out_best[i, 4] < out[j, 4]) {
	 		out_best[i, ] = out[j, ]
	 		}
	 	}
     }
   }
  
  # in out_best_best elimino eventuali soluzioni che hanno fit minore di soluzioni con minore nfp
  
  out_best_best = out_best[1, ,drop = F] # 
  if(tmp3>=2){ # caso limite in cui ho lo stesso nfp per tutte le righe in out -> out_best ha una sola riga
    for (i in 2:tmp3) {
      if (max(out_best_best[,4]) < out_best[i, 4]) {
        out_best_best = rbind(out_best_best, out_best[i, ])
      }
    }
  }
  tmp4 = nrow(out_best_best)
  
  # in out_best_best elimino eventuali soluzioni con punto sotto la linea
  if(tmp4>2){
    i = 2
    while (i < (tmp4 - 1)) {
	  # fit di tripla di soluzioni successive
      f1 = out_best_best[(i - 1), 4]	
      f2 = out_best_best[i, 4]
      f3 = out_best_best[(i + 1), 4]
      
	  # nfp di tripla di soluzioni successive
      fp1 = out_best_best[(i - 1), 6]	
      fp2 = out_best_best[i, 6]
      fp3 = out_best_best[(i + 1), 6]
      
	  # con LineCon verifico se sono sotto la linea
      if (LineCon(f1, f2, f3, fp1, fp2, fp3) == 0) {
	  	out_best_best = rbind(out_best_best[1:(i - 1), ], out_best_best[(i + 1):tmp4, ])
	  	tmp4 = nrow(out_best_best)
	  	i = 1
      }
      i = i + 1
    }
    tmp4 = nrow(out_best_best)
    st = cbind(vector(mode = "numeric", length = tmp4))
    for (i in 2:(tmp4 - 1)) {
	  fi = out_best_best[i, 4]
      fi_p = out_best_best[(i - 1), 4]
      fi_n = out_best_best[(i + 1), 4]
      fpi = out_best_best[i, 6]
      fpi_p = out_best_best[(i - 1), 6]	
      fpi_n = out_best_best[(i + 1), 6]
  
      if (st[i - 1] < Inf) {
	  	st[i] = ((fi - fi_p)/(fpi - fpi_p))/((fi_n - fi)/(fpi_n - fpi))
      }
      else {
	  	st[i] = Inf
      }
    }
  }
  if(tmp4==1){
    st = 0
  }
  if(tmp4==2){
    st = c()
    st = c(st, out_best_best[1,4]/out_best_best[1,6])
    st = c(st, out_best_best[2,4]/out_best_best[2,6])
  }
  out_st = cbind(out_best_best, st)
  
  colnames(out_st) <- c(colnames(out),"st")
  print(out_st)
  #print(out_st)
  ## imax = SUM(out_st)$valueMaxc[6] # sembra non andare
  st_max = out_st[which.max(out_st[,10]), ]
  # cat("Table: Goodness-of-fit values (%) and Scree test values st", fill = TRUE)
  # cat("of the solutions on the higher boundary of the convex hull", fill = TRUE)
  # tab = matrix(0, nrow = tmp4, ncol = 2)
  # for (i in 1:tmp4) {
	# tab[i, 1] = noquote(paste("     (", out_st[i, 1],                                 ") "))
	# if (i == 1) {
		# tab[i, 2] = noquote(paste("  ", round(out_st[i, 4], digits = 4), "      -"))
	# }
	# if (i == tmp4) {
		# tab[i, 2] = noquote(paste("  ", round(out_st[i, 4], digits = 4), "      -"))
	# }
	# else {
		# tab[i, 2] = noquote(paste("  ", round(out_st[i, 4], digits = 4), "    ", round(out_st[i, 6], digits = 2)))
	# }
	# if (i == imax) {
		# suggestion = noquote(paste("The hull heuristic indicates the selection of the CP model of complexity (", out_st[i, 1], ")"))
	# }
  # }
  # rownames(tab) = rep("CP", length = tmp4)
  # colnames(tab) = c("   Complexity", "   Fit (%)      st")
  # print(tab, quote = FALSE)
  # cat(suggestion, fill = TRUE)
  return(list(out_st = out_st, st_max = st_max))
}

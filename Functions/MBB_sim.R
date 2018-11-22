#Função para Moving Blocks Bootstrap

# function(numero de series, STL, Residuos, 
# tamanho do bloco, numero de dias, plotar graficos? (T ou F)))

MBB_function <- function(n_series,sinal,residuos,l,n_days,graf,label)
{
  n = length(residuos)
  #l = 48
  b = n - l + 1
  bloco = NULL
  for (i in 1:b) {
    bloco = cbind(bloco, residuos[i:(i + l - 1)])
  }
  blocos = as.matrix(bloco)
  simulacoes_boot_MBB = amostras_boot_MBB_ = matrix(NA, ncol = n_series, nrow = n)
  for (i in 1:n_series) {
    set.seed(i);selecao = sample(1:dim(bloco)[2], ceiling(n/l), replace = T)
    ruido.bs.prel = c(blocos[, selecao])
    ruido.bs = ruido.bs.prel[(length(ruido.bs.prel) - 
                                length(residuos) + 1):length(ruido.bs.prel)]
    simulacoes_boot_MBB[, i] = ruido.bs
    amostras_boot_MBB_[, i] = sinal[1:(24*n_days)] + simulacoes_boot_MBB[, i]
  }
  
  amostras_boot_MBB <- as.matrix(transpose(data.table(amostras_boot_MBB_)))
  amostras_boot_MBB_mean <- colMeans(amostras_boot_MBB)
  
  boot_distr <- array(dim = c(n_days*n_series,24))
  boot_mstl_distr <- array(dim = c(n_days*n_series,24))
  for (i in 1:24){
    for (s in 1:n_series){
      for (d in 1:n_days){
        
        boot_distr[(s-1)*n_days+d,i] <- amostras_boot_MBB[s,(d-1)*24+i] 
        boot_mstl_distr[(s-1)*n_days+d,i] <- amostras_boot_MBB[s,(d-1)*24+i] 
      }
    }
  }
  
  mstl_boot_mean <- colMeans(boot_mstl_distr)
  assign(paste(label,"_boot",sep=""),boot_mstl_distr, envir= .GlobalEnv)
  assign(paste(label,"_boot_mean",sep=""),mstl_boot_mean, envir= .GlobalEnv)
  assign(paste(label,"_boot_total",sep=""),amostras_boot_MBB_, envir= .GlobalEnv)
  
  if (graf == T){
    matplot(x = t(amostras_boot_MBB[,1:(24*n_days)]), type = "l", 
            lty = 1, col=alpha("gray",0.5), ylim = c(0,1),
            main = "Simulated Series: MSTL + MBB", 
            xlab = "Hour",
            ylab = "Capacity Factor")
    
    
    fan(amostras_boot_MBB[,1:(24*n_days)], ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
    lines(amostras_boot_MBB_mean[1:(24*n_days)], lw=2)
    #axis(side=1, at=c(0,1:24))
    
    #fanplot #2
    plot(mstl_boot_mean, type = "l",ylim = c(0,1), col="blue",ylab = "Capacity Factor",
         xlab = "Hour")
    fan(boot_mstl_distr,fan.col = sequential_hcl,ln=c(5, 25, 50, 75, 95))
    lines(mstl_boot_mean, type = "l",ylim = c(0,1), col="black")
    axis(side=1, at=c(0,1:24))
    
  }
  
}

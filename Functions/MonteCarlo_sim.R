# Função para simulação de MonteCarlo dos ruídos (dist. Normal)

remainder_mc <- function(n_series,sinal,mean_mc,sd_mc,day_ini,n_days,label,graf)
  #function(numero de series, série STL ou sinal, média dos ruidos, desvio padrao dos ruidos)
{
 
  len_days  <- n_days*24
   #n_series <- 1000
  # = 1 #Número de dias para considerar os cenários
  mc_mstl_series <- array(dim = c(n_series,len_days))
  monthly_mean_mstl_series <- array(dim = c(n_series,len_days))

  
  for (i in 1:n_series){
    set.seed(i^2);mc_mstl_series[i,] <- 
      sinal[((day_ini)*24-23):((day_ini)*24+(n_days-1)*24)] + 
      rnorm(n = len_days, mean_mc, sd_mc)
    mstl_series_matrix <- transpose(data.table(as.matrix(mc_mstl_series[i,])))
    monthly_mean_mstl_series[i,] <- colMeans(mstl_series_matrix)
  }
  

  mc_distr <- array(dim = c(n_days*n_series,24))
  #Distribuição de fatores de capacidade para cada hora do mês [31*n_series,24]
  mc_mstl_distr <- array(dim = c(n_days*n_series,24))
  for (i in 1:24){
    for (s in 1:n_series){
      for (d in 1:n_days){
        
        mc_distr[(s-1)*n_days+d,i] <- mc_mstl_series[s,(d-1)*24+i] 
        mc_mstl_distr[(s-1)*n_days+d,i] <- mc_mstl_series[s,(d-1)*24+i] 
      }
    }
  }
  
  monthly_expscen_mstl <- colMeans(mc_mstl_distr) #Valor esperado dos cenários para um dia típico
  
  assign(paste(label,"_montecarlo",sep=""),mc_mstl_distr, envir= .GlobalEnv)
  assign(paste(label,"_montecarlo_mean",sep=""),monthly_expscen_mstl,envir= .GlobalEnv)
  assign(paste(label,"_montecarlo_total",sep=""),mc_mstl_series, envir= .GlobalEnv)
  
  if(graf == TRUE){
    matplot(x = t(mc_mstl_series[,1:len_days]), type = "l", 
            lty = 1, col=alpha("gray",0.5), ylim = c(0,1),
            main = "Simulated Series: MSTL + Monte Carlo", 
            xlab = "Hour",
            ylab = "Capacity Factor")
    
    fan(mc_mstl_series, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
    lines(colMeans(mc_mstl_series),ylim = c(0,1), type = "l", col = "black", lw = 2)
    #axis(side=1, at=c(0,1:24))
    #Outras opções para fanplot
    plot(monthly_expscen_mstl, type = "l",ylim = c(0,1), col="blue",
         ylab = "Capacity Factor", xlab = "Hour")
    fan(mc_distr[,1:24],fan.col = sequential_hcl,ln=c(5, 25, 50, 75, 95),
        ln.col="red")
    axis(side=1, at=c(0,1:24))
    #fan(mc_distr,ln=c(5, 25, 50, 75, 95), style = "spaghetti",n.spag = 100)
    #fan(mc_mstl_distr, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
    }
}

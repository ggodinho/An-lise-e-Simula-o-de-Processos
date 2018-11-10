library(tseries)
library(TSA)
library(forecast)
library(MASS)
library(nortest)
library(data.table)
library(ggplot2)
library(xts)
library(fanplot)
library(boot)
library(fitdistrplus)
library(logspline)
library(colorspace)
# https://otexts.org/fpp2/forecasting-decomposition.html
###############################################################################
# FUNÁ√O PARA CONSTRU«√O DOS GR¡FICOS DE FORECAST ERRORS ######################
plotForecastErrors <- function(forecasterrors)
{
  # make a histogram of the forecast errors:
  mybinsize <- IQR(forecasterrors)/4
  mysd   <- sd(forecasterrors)
  mymin  <- min(forecasterrors) - mysd*5
  mymax  <- max(forecasterrors) + mysd*3
  # generate normally distributed data with mean 0 and standard deviation mysd
  mynorm <- rnorm(10000, mean=0, sd=mysd)
  mymin2 <- min(mynorm)
  mymax2 <- max(mynorm)
  if (mymin2 < mymin) { mymin <- mymin2 }
  if (mymax2 > mymax) { mymax <- mymax2 }
  # make a red histogram of the forecast errors, with the normally distributed data overlaid:
  mybins <- seq(mymin, mymax, mybinsize)
  hist(forecasterrors, col="red", freq=FALSE, breaks=mybins, xlab = "ResÌduos", ylab = "Densidade")
  # freq=FALSE ensures the area under the histogram = 1
  # generate normally distributed data with mean 0 and standard deviation mysd
  myhist <- hist(mynorm, plot=FALSE, breaks=mybins)
  # plot the normal curve as a blue line on top of the histogram of forecast errors:
  points(myhist$mids, myhist$density, type="l", col="blue", lwd=2)
}

######################################################
# Leitura dos dados #
######################################################
path <- getwd()
csv = read.csv2("./fcapacidade.csv",header=FALSE) #SÈrie HistÛrica de f. de capacidade
verif2018 <- as.matrix(read.csv2("./2018_verified.csv",header=F)) #Valores verificados em jan/2018
modelo_vigente <- as.matrix(read.csv2("./Modelo_Vigente.csv")) #SÈries usadas no modelo vigente

######################################################
# An·lise da DistribuiÁ„o da sÈrie #
######################################################
series <- ts(csv[,2], frequency = 24)
adf.test(series, alternative = "stationary") #SÈrie È estacion·ria, se rejeita hipÛtese nula

#https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# DistribuiÁıes candidatas:
descdist(csv[,2], discrete = FALSE) #Provavelmente Uniforme ou Normal
fit.normal <- fitdist(csv[,2],"norm")
plot(fit.normal)


######################################################
# Decompondo com MSTL #
######################################################
mstl_decomp <- mstl(series, iterate = 10, t.window = 730, s.window = 24)
plot(mstl_decomp)
plotForecastErrors(mstl_decomp[,4])
Box.test(mstl_decomp[,4], type = "Ljung-Box") #se p-valor <0.05, rejeito H0 (erros iid) com 5% de signific‚ncia

mstl_agg <- mstl_decomp[,2] + mstl_decomp[,3] #SÈrie STL (sem remainder)
mstl_agg_jan <- array(dim = c(31,24)) #SÈrie STL mÈdia de janeiro
for (d in 1:31){
  for (i in 1:24){
    mstl_agg_jan[d,i] <- mstl_agg[(d-1)*24+i] 
  }
}
mstl_agg_jan <- colMeans(mstl_agg_jan)
  
mstl_res <- mstl_decomp[,4] #Remainder
mean_mstl <- mean(mstl_decomp[,4]) #MÈdia ~ 0
sd_mstl <- sd(mstl_decomp[,4]) #Desvio padr„o

######################################################
# MSTL + MonteCarlo #
######################################################
remainder_mc <- function(n_series,sinal,mean_mc,sd_mc)
{
  #n_series <- 1000
  # = 1 #N˙mero de dias para considerar os cen·rios
  mc_mstl_series <- array(dim = c(n_series,24))
  monthly_mean_mstl_series <- array(dim = c(n_series,24))
  
  for (i in 1:n_series){
    set.seed(i^2);mc_mstl_series[i,] <- mstl_agg_jan + 
      rnorm(n = 24, mean_mc, sd_mc)
    mstl_series_matrix <- transpose(data.table(matrix(mc_mstl_series[i,],24,31)))
    monthly_mean_mstl_series[i,] <- colMeans(mstl_series_matrix)
    lines(mc_mstl_series[i,1:24],ylim = c(0,1), type = "l", 
          col=alpha("gray",0.5))
  }
  
  fan(mc_mstl_series[,1:24], ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
  lines(colMeans(mc_mstl_series[,1:24]),ylim = c(0,1), type = "l", col = "black", lw = 2)
  axis(side=1, at=c(0,1:24))
  monthly_expscen_mstl <- colMeans(monthly_mean_mstl_series) #Valor esperado dos cen·rios para um dia tÌpico
  
  assign("montecarlo_series",mc_mstl_series, envir= .GlobalEnv)
  assign("montecarlo_series_mean",monthly_expscen_mstl,envir= .GlobalEnv)
  
  #Outras opÁıes para fanplot
  plot(monthly_expscen_mstl, type = "l",ylim = c(0,1), col="blue",
       ylab = "Capacity Factor", xlab = "Hour")
  fan(mc_mstl_series,fan.col = sequential_hcl,ln=c(5, 25, 50, 75, 95))
  axis(side=1, at=c(0,1:24))
  #fan(mc_distr,ln=c(5, 25, 50, 75, 95), style = "spaghetti",n.spag = 100)
  #fan(mc_mstl_distr, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
}

remainder_mc(1000,mstl_agg_jan,0,sd_mstl) #Monta 2 gr·ficos

######################################################
# MSTL + Bootstrap #
######################################################

MBB_function <- function(n_series,sinal,residuos,l)
{
  n = length(residuos)
  #l = 48
  b = n - l + 1
  bloco = NULL
  for (i in 1:b) {
    bloco = cbind(bloco, residuos[i:(i + l - 1)])
  }
  blocos = as.matrix(bloco)
  simulacoes_boot_MBB = amostras_boot_MBB = matrix(NA, ncol = n_series, nrow = n)
  for (i in 1:n_series) {
    set.seed(i);selecao = sample(1:dim(bloco)[2], ceiling(n/l), replace = T)
    ruido.bs.prel = c(blocos[, selecao])
    ruido.bs = ruido.bs.prel[(length(ruido.bs.prel) - length(residuos) + 1):length(ruido.bs.prel)]
    simulacoes_boot_MBB[, i] = ruido.bs
    amostras_boot_MBB[, i] = sinal + simulacoes_boot_MBB[, i]
  }
  
  amostras_boot_MBB <- as.matrix(transpose(data.table(amostras_boot_MBB[1:24,])))
  mstl_boot_mean <- colMeans(amostras_boot_MBB)
  plot(mstl_boot_mean[1:24],type = "l",ylim = c(0,1),ylab = "Capacity Factor",
       xlab = "Hour")
  for (i in 1:n_series){
    lines(amostras_boot_MBB[i,1:24],col=alpha("gray",0.5))
  }

  fan(amostras_boot_MBB, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
  lines(mstl_boot_mean, lw=2)
  axis(side=1, at=c(0,1:24))
  
  #fanplot #2
  plot(mstl_boot_mean, type = "l",ylim = c(0,1), col="blue",ylab = "Capacity Factor",
       xlab = "Hour")
  fan(amostras_boot_MBB,fan.col = sequential_hcl,ln=c(5, 25, 50, 75, 95))
  lines(mstl_boot_mean, type = "l",ylim = c(0,1), col="black")
  axis(side=1, at=c(0,1:24))
  
  assign("amostras_boot_MBB",amostras_boot_MBB, envir= .GlobalEnv)
  assign("amostras_boot_MBB_mean",mstl_boot_mean, envir= .GlobalEnv)
}

MBB_function(1000,mstl_agg_jan,mstl_res,l = 48) #Monta 2 gr·ficos

######################################################
# Comparando com Modelo Vigente #
######################################################
prob_curmodel = modelo_vigente[,26]

exp_curmodel = (modelo_vigente[1,2:25]*prob_curmodel[1]+
                  modelo_vigente[2,2:25]*prob_curmodel[2]+
                  modelo_vigente[3,2:25]*prob_curmodel[3])

# Current Model x MSTL+MonteCarlo
# Comparando intervalos de confianÁa x MSTL+MonteCarlo
plot(montecarlo_series_mean, type = "l",ylim = c(0,1), col="blue",ylab = "Capacity Factor",
     xlab = "Hour",main = "Current Model x MSTL+MonteCarlo")
fan(montecarlo_series,fan.col = sequential_hcl,ln=c(5, 95),lw=2)
lines(modelo_vigente[1,2:25],col="dark red",lw=3)
lines(modelo_vigente[2,2:25],col="red",lw=3)
lines(modelo_vigente[3,2:25],col="orange",lw=3)
axis(side=1, at=c(0,1:24))

#Comparando Valor Esperado dos Cen·rios
plot(exp_curmodel, type = "l", ylim = c(0,1), col="red",ylab = "Capacity Factor",
     xlab = "Hour", main = "Expected Value: Current Model x MSTL+MonteCarlo", lwd=2)
lines(montecarlo_series_mean,col="blue",lw=2,lty=2)
axis(side=1, at=c(0,1:24))
legend("bottomright",legend = c("Current Model","MSTL+MonteCarlo"),
       col=c("red","blue"), lwd = c(2,2), lty= c(1,2), bty="n")


# Current Model x MSTL+Bootstrap
# Comparando intervalos de confianÁa x MSTL+Bootstrap
plot(amostras_boot_MBB_mean, type = "l",ylim = c(0,1), col="blue",ylab = "Capacity Factor",
     xlab = "Hour",main = "Current Model x MSTL+Bootstrap")
fan(amostras_boot_MBB,fan.col = sequential_hcl,ln=c(5, 95),lw=2)
lines(modelo_vigente[1,2:25],col="dark red",lw=3)
lines(modelo_vigente[2,2:25],col="red",lw=3)
lines(modelo_vigente[3,2:25],col="orange",lw=3)
axis(side=1, at=c(0,1:24))

#Comparando Valor Esperado dos Cen·rios
plot(exp_curmodel, type = "l", ylim = c(0,1), col="red",ylab = "Capacity Factor",
     xlab = "Hour", main = "Expected Value: Current Model x MSTL+Bootstrap", lwd=2)
lines(amostras_boot_MBB_mean,col="blue",lw=2,lty=2)
axis(side=1, at=c(0,1:24))
legend("bottomright",legend = c("Current Model","MSTL+Bootstrap"),
       col=c("red","blue"), lwd = c(2,2), lty= c(1,2), bty="n")

# Comparando a densidade de distribuiÁ„o das sÈries
d1=density(verif2018)
d2=density(amostras_boot_MBB)           
d3=density(modelo_vigente[,2:25])
           
plot(d1,ylim=c(0,6), lwd = 1, xlab = "Capacity Factor",xlim=c(0,1),
     main="Density")
polygon(d1,col=alpha("gray",0.5))
lines(d2, col = "blue",lty=2,lwd=2)
lines(d3, col = "red",lwd=2)
legend("topright", legend = c("Time Series","MSTL+Bootstrap","Current Model"),
       bty="n", fill=c(alpha("gray",0.5), NA,NA),xpd=T,lty=c(NA,2,1),
       border = c("light blue", NA,NA),horiz=F,col=c(0,"blue","red"),
       pch = c(0,-1,-1))

######################################################
# Comparando com a sÈrie verificada out-of-sample #
######################################################

verif2018_m <- array(dim = c(31,24)) 
for (d in 1:31){
  for (i in 1:24){
    verif2018_m[d,i] <- verif2018[(d-1)*24+i] 
  }
}
verif2018_mean <- colMeans(verif2018_m)
plot(verif2018,type="l",ylim=c(0,1))

plot(verif2018_m[1,], type = "l",ylim = c(0,1),col="gray",xlab="Hour",
     ylab="Capacity Factor",main = "Out-of-sample data x Models Used")
for(i in 2:31){
  lines(verif2018_m[i,],col = "gray")
}
lines(verif2018_mean,type = "l", ylim = c(0,1),lwd=2, col = alpha("black",0.5))
lines(exp_curmodel,type = "l", ylim = c(0,1),lwd=2, col = "red")
lines(montecarlo_series_mean,type = "l", ylim = c(0,1),lwd=2, col = "blue")
lines(amostras_boot_MBB_mean,type = "l", ylim = c(0,1),lwd=2, col = "green")
legend("topright", legend = c("Out-of-sample data","Current Model",
                              "MSTL+MonteCarlo","MSTL+Bootstrap"),
       bty="n", xpd=T,lty=c(1,1,1,1), horiz=F,
       col=c(alpha("black",0.5),"red","blue","green"),lwd = c(2,2,2,2))

######################################################
# Calculando MAPE
######################################################

#MAPE total
mape_curmodel = mean(abs((exp_curmodel-verif2018_mean)/(verif2018_mean)))*100
mape_mstl_montecarlo = mean(abs((montecarlo_series_mean-verif2018_mean)/(verif2018_mean)))*100
mape_mstl_boot = mean(abs((amostras_boot_MBB_mean-verif2018_mean)/(verif2018_mean)))*100

#MAPE hor·rio
windowmape <- function(model_mean,verified_mean,window)
{
  aux = NULL
  for(i in 1:window){
    aux[i] = accuracy(model_mean[i],verified_mean[i])[5]
  }
  assign("monthlymape",aux, envir= .GlobalEnv)
}

hmape_curmodel <- windowmape(exp_curmodel,verif2018_mean,24)
hmape_mstl_montecarlo <- windowmape(montecarlo_series_mean,verif2018_mean,24)
hmape_mstl_boot <- windowmape(amostras_boot_MBB_mean,verif2018_mean,24)

plot(hmape_curmodel,type="l",lwd=2,xlab="Hour", ylab="MAPE (%)",ylim = c(0,40), 
     main = "Hourly MAPE",col = "red")
lines(hmape_mstl_montecarlo,type = "l",lwd=2, col = "blue")
lines(hmape_mstl_boot,type = "l", lwd=2, col = "green")
legend("topright", legend = c("Current Model","MSTL+MonteCarlo","MSTL+Bootstrap"),
       bty="n", xpd=T,lty=c(1,1,1), horiz=F,
       col=c("red","blue","green"),lwd = c(2,2,2))

######################################################
# Testando ETS #
######################################################
# ets_decomp <- ets(series)
# plot(ets_decomp)
# set.seed(0)
# 
# ets_series <- array(dim = c(n_series,744))
# for (k in 1:n_series){
#   ets_series[k,] <- simulate(ets_decomp, nsim = 744, seed = k)
# }
# 
# plot(ets_series[1,], ylim = c(0,1),col = "gray", type = "line")
# for (k in 2:n_series){
#   lines(ets_series[k,], ylim = c(0,1),col = "gray", type = "line")
# }
# 
# ets_series_mean <- colMeans(ets_series)
# lines(ets_series_mean, ylim = c(0,1),col = "red", type = "line",lwd=2)
# 
# 
# ######################################################
# # Ajustando um modelo SARIMA #
# ######################################################
# 
# fit_arima <- auto.arima(series)
# summary(fit_arima)
# tsdisplay(residuals(fit_arima))
# 
# fit_arima2 <- Arima(series,order=c(4,1,1),seasonal=list(order=c(2,0,1),period=24))
# summary(fit_arima2)
# tsdisplay(residuals(fit_arima2))
# 
# forecast <- forecast(seasonaldecomp, h = 8760)
# plot(forecast)
# 
# ##
# ##save as png
# ##
# # dev.copy(png, file = "svplots.png", width=10, height=50, units="in", res=400)
# # dev.off()

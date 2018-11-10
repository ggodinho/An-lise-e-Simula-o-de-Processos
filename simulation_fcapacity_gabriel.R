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
#Leitura dos dados
path <- getwd()
csv = read.csv2("./fcapacidade.csv",header=FALSE)

## An·lise da DistribuiÁ„o da sÈrie ##################
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
mstl_res <- mstl_decomp[,4] #Remainder
mean_mstl <- mean(mstl_decomp[,4]) #MÈdia ~ 0
sd_mstl <- sd(mstl_decomp[,4]) #Desvio padr„o

#Simulando Monte Carlo
remainder_mc <- function(n_series,sinal,n_days,mean_mc,sd_mc)
{
  #n_series <- 100
  # = 1 #N˙mero de dias para considerar os cen·rios
  mc_mstl_series <- array(dim = c(n_series,n_days*24))
  monthly_mean_mstl_series <- array(dim = c(n_series,24))
  
  plot(mstl_agg[1:(24*n_days)],ylim = c(0,1), type = "l", col = "red")
  for (i in 1:n_series){
    set.seed(i^2);mc_mstl_series[i,] <- mstl_agg[1:(24*n_days)] + 
      rnorm(n = (n_days*24), mean_mc, sd_mc)
    mstl_series_matrix <- transpose(data.table(matrix(mc_mstl_series[i,],24,31)))
    monthly_mean_mstl_series[i,] <- colMeans(mstl_series_matrix)
    lines(mc_mstl_series[i,1:(24*n_days)],ylim = c(0,1), type = "l", 
          col=alpha("gray",0.5))
  }
  
  fan(mc_mstl_series[,1:(24*n_days)], ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
  lines(colMeans(mc_mstl_series[,1:(24*n_days)]),ylim = c(0,1), type = "l", col = "black", lw = 2)
  axis(side=1, at=c(0,1:(n_days*24)))
  
  monthly_expscen_mstl <- colMeans(monthly_mean_mstl_series) #Valor esperado dos cen·rios para um dia tÌpico
  
  #Outras opÁıes para fanplot
  #plot(monthly_expscen_mstl, type = "l",ylim = c(0,1), col="blue")
  #fan(mc_distr,fan.col = sequential_hcl,ln=c(5, 25, 50, 75, 95))
  #fan(mc_distr,ln=c(5, 25, 50, 75, 95), style = "spaghetti",n.spag = 100)
  #fan(mc_mstl_distr, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
}

remainder_mc(500,mstl_agg,1,0,sd_mstl)

# Fazer bootstrap dos ruidos para MSTL

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
  plot(mstl_boot_mean[1:24],type = "l",ylim = c(0,1))
  for (i in 1:n_series){
    lines(amostras_boot_MBB[i,1:24],col=alpha("gray",0.5))
  }
  fan(amostras_boot_MBB, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
  lines(mstl_boot_mean, lw=2)
}

MBB_function(500,mstl_agg,mstl_res,l = 48)



# mstl_res_boot <- as.matrix(transpose(data.table(tsboot(mstl_res, 
#                                                   n.sim=1000, l=48))))
# mstl_res_boot <- mstl_res_boot[,1:24]
# 
# mstl_boot <- mstl_agg[1:24] + mstl_res_boot
# mstl_boot_mean <- colMeans(mstl_boot)

plot(mstl_boot_mean,type = "l",ylim = c(0,1))
for (i in 1:n_series){
  lines(mstl_boot[i,],col=alpha("gray",0.5))
}
fan(mstl_boot[,1:(24*n_days)], ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
lines(mstl_boot_mean, lw=2)

plot(mstl_agg[1:24],ylim=c(0,1))

#Comparar mÈdias e QQPLOT da distribuiÁ„o verificada em 2018 e distribuiÁ„o simulada 

######################################################
# Testando ETS #
######################################################
ets_decomp <- ets(series)
plot(ets_decomp)
set.seed(0)

ets_series <- array(dim = c(n_series,744))
for (k in 1:n_series){
  ets_series[k,] <- simulate(ets_decomp, nsim = 744, seed = k)
}

plot(ets_series[1,], ylim = c(0,1),col = "gray", type = "line")
for (k in 2:n_series){
  lines(ets_series[k,], ylim = c(0,1),col = "gray", type = "line")
}

ets_series_mean <- colMeans(ets_series)
lines(ets_series_mean, ylim = c(0,1),col = "red", type = "line",lwd=2)


######################################################
# Ajustando um modelo SARIMA #
######################################################

fit_arima <- auto.arima(series)
summary(fit_arima)
tsdisplay(residuals(fit_arima))

fit_arima2 <- Arima(series,order=c(4,1,1),seasonal=list(order=c(2,0,1),period=24))
summary(fit_arima2)
tsdisplay(residuals(fit_arima2))

forecast <- forecast(seasonaldecomp, h = 8760)
plot(forecast)

##
##save as png
##
# dev.copy(png, file = "svplots.png", width=10, height=50, units="in", res=400)
# dev.off()

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

mean_mstl <- mean(mstl_decomp[,4]) #MÈdia ~ 0
sd_mstl <- sd(mstl_decomp[,4]) #Desvio padr„o
mstl_agg <- mstl_decomp[,2] + mstl_decomp[,3] #SÈrie STL (sem remainder)

#Simulando Monte Carlo
n_series <- 100
mc_mstl_series <- array(dim = c(n_series,744))
monthly_mean_mstl_series <- array(dim = c(n_series,24))
n_days = 1 #Quantos dias plotar

plot(mstl_agg[1:(24*n_days)],ylim = c(0,1), type = "l", col = "red")
for (i in 1:n_series){
  set.seed(i^2);mc_mstl_series[i,] <- mstl_agg[1:(24*n_days)] + 
    rnorm(n = 744, mean = 0, sd = sd_mstl)
  mstl_series_matrix <- transpose(data.table(matrix(mc_mstl_series[i,],24,31)))
  monthly_mean_mstl_series[i,] <- colMeans(mstl_series_matrix)
  lines(mc_mstl_series[i,1:(24*n_days)],ylim = c(0,1), type = "l", 
        col=alpha("gray",0.5))
}

fan(mc_mstl_series[,1:(24*n_days)], ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
lines(colMeans(mc_mstl_series[,1:(24*n_days)]),ylim = c(0,1), type = "l", col = "black", lw = 2)
axis(side=1, at=c(0,1:24))


monthly_expscen_mstl <- colMeans(monthly_mean_mstl_series) #Valor esperado dos cen·rios para um dia tÌpico

#DistribuiÁ„o de fatores de capacidade para cada hora do mÍs [31*n_series,24]
mc_mstl_distr <- array(dim = c(31*n_series,24))
for (i in 1:24){
  for (s in 1:n_series){
    for (d in 1:31){
      
      mc_mstl_distr[(s-1)*31+d,i] <- mc_mstl_series[s,(d-1)*24+i] 
    }
  }
}
fan(mc_mstl_distr, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")
lines(monthly_expscen_mstl, type = "l",ylim = c(0,1), col="black",lw = 3)
axis(side=1, at=c(0,1:24))

# Fan plot com a distribuiÁ„o em um dia tÌpico de janeiro
plot(monthly_expscen_mstl, type = "l",ylim = c(0,1), col="blue")
for (i in 1:(31*n_series)){
  lines(mc_mstl_distr[i,],ylim = c(0,1), type = "l", col=alpha("blue",0.02))
}


#Outras opÁıes para fanplot
#plot(monthly_expscen_mstl, type = "l",ylim = c(0,1), col="blue")
#fan(mc_distr,fan.col = sequential_hcl,ln=c(5, 25, 50, 75, 95))
#fan(mc_distr,ln=c(5, 25, 50, 75, 95), style = "spaghetti",n.spag = 100)
#fan(mc_mstl_distr, ln=c(5, 25, 50, 75, 95), alpha=0,ln.col="red")

# Fazer bootstrap dos ruidos para MSTL

######################################################
# MSTL + Bootstrap #
######################################################

mstl_boot <- array(dim = c(n_series,length(series)))

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

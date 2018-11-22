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
source("./Functions/plotForecastErrors.R")
source("./Functions/MonteCarlo_sim.R")
source("./Functions/MBB_sim2.R")
source("./Functions/Benchmark_comparison.R")
source("./Functions/observedvalues_comparison.R")

######################################################
# Leitura dos dados #
######################################################
path <- getwd()
csv = read.csv2("./Data/fcapacidade.csv",header=FALSE) #Série Histórica de f. de capacidade
verif2018 <- as.matrix(read.csv2("./Data/2018_verified.csv",header=F)) #Valores verificados em jan/2018
modelo_vigente <- as.matrix(read.csv2("./Data/Modelo_Vigente.csv")) #Séries usadas no modelo vigente

######################################################
# Análise da Distribuição da série #
######################################################
series <- msts(csv[,2], seasonal.periods = c(24,8760))
adf.test(series, alternative = "explosive",k=730) #H0 (estacionaria) não rejeitada com 5% de significância. p-valor = 0,099

#https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# Distribuições candidatas:
descdist(csv[,2], discrete = FALSE) #Provavelmente Uniforme ou Normal
fit.normal <- fitdist(csv[,2],"norm")
plot(fit.normal)
ks.test(x = series, y = "pnorm")


######################################################
# Decompondo com MSTL #
######################################################
mstl_decomp <- mstl(series, iterate = 100, t.window = 730, s.window = 365)
plot(mstl_decomp, main = "MSTL Series Decomposition")
plot(mstl_decomp[1:24,3], type = "l")

mstl_agg <- mstl_decomp[,2] + mstl_decomp[,3] #Série STL (sem remainder)
mstl_agg_jan <- array(dim = c(31,24)) #Série STL média de janeiro
for (d in 1:31){
  for (i in 1:24){
    mstl_agg_jan[d,i] <- mstl_agg[(d-1)*24+i] 
  }
}
mstl_agg_jan <- colMeans(mstl_agg_jan)

mstl_res <- mstl_decomp[,4] #Remainder
mean_mstl <- mean(mstl_decomp[,4]) #Média ~ 0
sd_mstl <- sd(mstl_decomp[,4]) #Desvio padrão

descdist(as.vector(mstl_res), discrete = FALSE) #Provavelmente Uniforme ou Normal
fit.res <- fitdist(as.vector(mstl_res),"norm")
plot(fit.res)
plotForecastErrors(mstl_decomp[,4])
Box.test(mstl_decomp[,4], type = "Ljung-Box") #se p-valor <0.05, rejeito H0 (erros iid) com 5% de significância
tsdisplay(mstl_decomp[,4])
ks.test(x = mstl_decomp[,4], y = "pnorm")


######################################################
# MSTL + MonteCarlo #
######################################################

remainder_mc(1000,mstl_agg,mean_mstl,sd_mstl,day_ini = 1,n_days = 31,"mstl",
             graf = T) #Monta 2 gráficos

######################################################
# MSTL + Bootstrap #
######################################################

MBB_function(1000,mstl_agg,mstl_res,l = 48,n_days = 31,graf = T,
             label = "mstl") #Monta 2 gráficos

######################################################
# ARIMA #
######################################################
#AutoArima
fit_arima <- auto.arima(series)
summary(fit_arima)
tsdisplay(residuals(fit_arima),lag.max = 100)


#Ajustado - sobrefixacao
fit_arima2 <- Arima(series,order=c(4,1,1),seasonal=list(order=c(2,0,1),period=24))
summary(fit_arima2)
tsdisplay(residuals(fit_arima2),lag.max=100)
plotForecastErrors(residuals(fit_arima2))

#simArima = arima.sim(model = fit_arima, n = 100, rand.gen = rnorm) #arima.sim não simula modelos 

cenarioHorizonte = 24
nCen = 1000

# Cadeia de Markov --------------------------------------------------------
simArima <- array(dim = c(nCen,24))

for (i in 1:nCen){
  simArima[i,] <- simulate(fit_arima2, nsim = 24, seed = i, future = TRUE)
}

matplot(x = t(simArima), type = "l", lty = 1, col = "grey", ylim = c(0,1),
        main = "Simulated Series: ARIMA + Monte Carlo", 
        xlab = "Hour",
        ylab = "Capacity Factor")
lines(apply(simArima, 2, mean), lwd = 2, col = "red")
axis(side=1, at=c(0,1:24))


# ARIMA + Bootstrap ---------------------------------------------------------------
simArima_bt <- array(dim = c(nCen,24))
for (i in 1:nCen){
  simArima_bt[i,] <- simulate(fit_arima2, nsim = 24, bootstrap = T, 
                              innov = NULL, seed = i)
}

matplot(x = t(simArima_bt), type = "l", lty = 1, col = "grey", ylim = c(0,1),
        main = "Simulated Series: ARIMA + Bootstrap", 
        xlab = "Hour",
        ylab = "Capacity Factor")
lines(apply(simArima_bt, 2, mean), lwd = 2, col = "red")
axis(side=1, at=c(0,1:24))


######################################################
# Comparando com Modelo Vigente #
######################################################
prob_modelo_vigente = modelo_vigente[,26]

#Valor esperado Modelo Vigente
modelo_vigente_exp = (modelo_vigente[1,2:25]*prob_modelo_vigente[1]+
                  modelo_vigente[2,2:25]*prob_modelo_vigente[2]+
                  modelo_vigente[3,2:25]*prob_modelo_vigente[3])

#Comparando MSTL+MonteCarlo
benchmark_comp(mstl_montecarlo,mstl_montecarlo_mean,modelo_vigente,
               modelo_vigente_exp,"MSTL+MonteCarlo")

#Comparando MSTL+Bootstrap
benchmark_comp(mstl_boot,mstl_boot_mean,modelo_vigente,
               modelo_vigente_exp,"MSTL+Bootstrap")

#Comparando ARIMA+MonteCarlo
benchmark_comp(simArima,apply(simArima, 2, mean),modelo_vigente,
               modelo_vigente_exp,"ARIMA+MonteCarlo")

#Comparando ARIMA+Bootstrap
benchmark_comp(simArima_bt,apply(simArima_bt, 2, mean),modelo_vigente,
               modelo_vigente_exp,"ARIMA+Bootstrap")

######################################################
# Comparando com a série verificada out-of-sample #
######################################################

verif2018_m <- array(dim = c(31,24)) 
for (d in 1:31){
  for (i in 1:24){
    verif2018_m[d,i] <- verif2018[(d-1)*24+i] 
  }
}
verif2018_mean <- colMeans(verif2018_m)
plot(verif2018,type="l",ylim=c(0,1),main = "Out-of-sample data - January/2018",
     ylab = "Capacity Factor",xlab = "Hour")

observedvalues_comp(verif2018_m,modelo_vigente,"Current Model")
observedvalues_comp(verif2018_m,mstl_montecarlo,"MSTL+MonteCarlo")
observedvalues_comp(verif2018_m,mstl_boot,"MSTL+Bootstrap")
observedvalues_comp(verif2018_m,simArima,"ARIMA+MonteCarlo")
observedvalues_comp(verif2018_m,simArima_bt,"ARIMA+Bootstrap")

#Comparação Todos Modelos
#plot(verif2018_m[1,], type = "l",ylim = c(0,1),col="gray",xlab="Hour",
     #ylab="Capacity Factor",main = "Observed Values 2018 x Simulation Models")
#for(i in 2:31){
#  lines(verif2018_m[i,],col = "gray")
#}

verif_2017_m <- array(dim = c(31,24))

for (i in 1:24){
  for (d in 1:31){
      
    verif_2017_m[d,i] <- series[(d-1)*24+i] 
  }
}

plot(apply(verif_2017_m, 2, mean),type = "l", ylim = c(0,1),lwd=2, lty = 2,
     col = "black", ylab="Capacity Factor",xlab="Hour",
     main = "Observed Values 2017 x Simulation Models")
lines(modelo_vigente_exp,type = "l", ylim = c(0,1),lwd=2, col = "red")
lines(mstl_montecarlo_mean,type = "l", ylim = c(0,1),lwd=2, col = "blue")
lines(mstl_boot_mean,type = "l", ylim = c(0,1),lwd=2, col = "green")
lines(apply(simArima, 2, mean),type = "l", ylim = c(0,1),lwd=2, col = "purple")
lines(apply(simArima_bt, 2, mean),type = "l", ylim = c(0,1),lwd=2, col = "pink")
legend("topright", legend = c("Observed Values 2017","Current Model",
                              "MSTL+MonteCarlo","MSTL+Bootstrap",
                              "Arima+Montecarlo","Arima+Bootstrap"),
       bty="n", xpd=T,lty=c(2,1,1,1,1,1), horiz=F,
       col=c("black","red","blue","green","purple","pink"),lwd = c(2,2,2,2,2,2))


######################################################
# Calculando MAPE
######################################################

#MAPE total
MAPE <- function(model_mean,observed_mean,label){
  mape = mean(abs((model_mean-observed_mean)/(observed_mean)))*100
  name = paste("mape_",label,sep="")
  assign(name,mape,envir = .GlobalEnv)
}

MAPE(modelo_vigente_exp,verif2018_mean,"curmodel")
MAPE(mstl_montecarlo_mean,verif2018_mean,"mstl_montecarlo")
MAPE(mstl_boot_mean,verif2018_mean,"mstl_boot")
MAPE(apply(simArima, 2, mean),verif2018_mean,"arima_montecarlo")
MAPE(apply(simArima_bt, 2, mean),verif2018_mean,"arima_boot")

#MAPE horário
windowmape <- function(model_mean,verified_mean,window)
{
  aux = NULL
  for(i in 1:window){
    aux[i] = accuracy(model_mean[i],verified_mean[i])[5]
  }
  assign("monthlymape",aux, envir= .GlobalEnv)
}

hmape_curmodel <- windowmape(modelo_vigente_exp,verif2018_mean,24)
hmape_mstl_montecarlo <- windowmape(mstl_montecarlo_mean,verif2018_mean,24)
hmape_mstl_boot <- windowmape(mstl_boot_mean,verif2018_mean,24)
hmape_arima_mc <- windowmape(apply(simArima, 2, mean),verif2018_mean,24)
hmape_arima_boot <- windowmape(apply(simArima_bt, 2, mean),verif2018_mean,24)


plot(hmape_curmodel,type="l",lwd=2,xlab="Hour", ylab="MAPE (%)",ylim = c(0,20), 
     main = "Hourly MAPE",col = "red")
lines(hmape_mstl_montecarlo,type = "l",lwd=2, col = "blue")
lines(hmape_mstl_boot,type = "l", lwd=2, col = "cyan", lty=2)
lines(hmape_arima_mc,type = "l",lwd=2, col = "dark green")
lines(hmape_arima_boot,type = "l", lwd=2, col = "green", lty=2)
axis(side=1, at=c(0,1:24))
legend("topright", legend = c("Current Model","MSTL+MonteCarlo","MSTL+Bootstrap",
                              "ARIMA+MonteCarlo","ARIMA+Bootstrap"),
       bty="n", xpd=T,lty=c(1,1,2,1,2), horiz=F,
       col=c("red","blue","cyan","dark green","green"),lwd = c(2,2,2))

######################################################
# Decompondo via TBATS #
######################################################
# Ficou pior do que a decomposição MSTL

##tbats_decomp <- tbats(series)
#plot(tbats_decomp) #decomposição da série
#tbats_components <- tbats.components(tbats_decomp)
#
#plot(tbats_decomp$errors) #resíduos
#sd_tbats <- sd(tbats_decomp$errors)
#tsdisplay(tbats_decomp$errors, lag.max = 50)
#
#Box.test(tbats_decomp$errors, type = "Ljung-Box") #se p-valor <0.05, rejeito H0 (erros iid) com 5% de significância
#tbats_dist <- fitdist(as.vector(tbats_decomp$errors),"norm")
#plot(tbats_dist)
#
#tbats_signal <- tbats_components[,2] + tbats_components[,3]+tbats_components[,4]
#
#remainder_mc(100,tbats_signal,0,sd_tbats,"TBATS")
#MBB_function(100,tbats_signal,tbats_decomp$errors,48,31,T,"TBATS")
#benchmark_comp(TBATS_boot,TBATS_boot_mean,
#               modelo_vigente,exp_curmodel,"TBATS+Bootstrap")
#MAPE(TBATS_boot_mean,verif2018_mean,"tbats_boot")

######################################################
# Ajustando ARIMAs em cima da série MSTL+Bootstrap #
######################################################

#cenarioHorizonte = 24
#nCen = 50
#
#sim_mstl_Arima <- array(dim = c(nCen,24))
#
#for (i in 1:nCen){
#  fit_mstl_Arima <- auto.arima(mstl_365_montecarlo_total[i,])
#  sim_mstl_Arima[i,] <- simulate(fit_mstl_Arima, nsim = 24, seed = i, 
#                                 future = TRUE)
#  #fc_mstl_Arima <- forecast(fit_mstl_Arima, h=24)
#  #sim_mstl_Arima[i,] <- fc_mstl_Arima$mean[1:24]
#}
#
#matplot(x = t(sim_mstl_Arima), type = "l", lty = 1, col = "grey", ylim = c(0,1),
#        main = "MSTL+ARIMA Simulations", 
#        xlab = "Horizonte horário",
#        ylab = "Fator de Carga para o NE")
#lines(apply(sim_mstl_Arima, 2, mean), lwd = 2, col = "red")


# 
# ##
# ##save as png
# ##
# # dev.copy(png, file = "svplots.png", width=10, height=50, units="in", res=400)
# # dev.off()

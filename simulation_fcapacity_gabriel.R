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
source("./Functions/MBB_sim.R")
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
adf.test(series, alternative = "explosive",k=744) #H0 (estacionaria) não rejeitada com 5% de significânica. p-valor = 0,099

#https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
# Distribuições candidatas:
descdist(csv[,2], discrete = FALSE) #Provavelmente Uniforme ou Normal
fit.normal <- fitdist(csv[,2],"norm")
plot(fit.normal)


######################################################
# Decompondo com MSTL #
######################################################
mstl_decomp <- mstl(series, iterate = 10, t.window = 730)
plot(mstl_decomp)

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


######################################################
# MSTL + MonteCarlo #
######################################################

remainder_mc(100,mstl_agg,0,sd_mstl) #Monta 2 gráficos

######################################################
# MSTL + Bootstrap #
######################################################

MBB_function(100,mstl_agg,mstl_res,l = 48,31,T) #Monta 2 gráficos

######################################################
# Comparando com Modelo Vigente #
######################################################
prob_curmodel = modelo_vigente[,26]

exp_curmodel = (modelo_vigente[1,2:25]*prob_curmodel[1]+
                  modelo_vigente[2,2:25]*prob_curmodel[2]+
                  modelo_vigente[3,2:25]*prob_curmodel[3])

#Comparando MSTL+MonteCarlo
benchmark_comp(montecarlo_series,montecarlo_series_mean,modelo_vigente,
               exp_curmodel,"MSTL+MonteCarlo")

#Comparando MSTL+Bootstrap
benchmark_comp(amostras_boot_MBB,amostras_boot_MBB_mean,modelo_vigente,
               exp_curmodel,"MSTL+Bootstrap")


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
plot(verif2018,type="l",ylim=c(0,1))

observedvalues_comp(verif2018_m,modelo_vigente,"Current Model")
observedvalues_comp(verif2018_m,montecarlo_series,"MSTL+MonteCarlo")
observedvalues_comp(verif2018_m,amostras_boot_MBB,"MSTL+Bootstrap")

#Comparação Todos Modelos
plot(verif2018_m[1,], type = "l",ylim = c(0,1),col="gray",xlab="Hour",
     ylab="Capacity Factor",main = "Observed Values 2018 x Simulation Models")
for(i in 2:31){
  lines(verif2018_m[i,],col = "gray")
}
lines(verif2018_mean,type = "l", ylim = c(0,1),lwd=2, col = alpha("black",0.5))
lines(exp_curmodel,type = "l", ylim = c(0,1),lwd=2, col = "red")
lines(montecarlo_series_mean,type = "l", ylim = c(0,1),lwd=2, col = "blue")
lines(amostras_boot_MBB_mean,type = "l", ylim = c(0,1),lwd=2, col = "green")
legend("topright", legend = c("Observed Values 2018","Current Model",
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

#MAPE horário
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
axis(side=1, at=c(0,1:24))
legend("topright", legend = c("Current Model","MSTL+MonteCarlo","MSTL+Bootstrap"),
       bty="n", xpd=T,lty=c(1,1,1), horiz=F,
       col=c("red","blue","green"),lwd = c(2,2,2))

######################################################
# Decompondo via TBATS #
######################################################
# Ficou pior do que a decomposição MSTL

#tbats_decomp <- tbats(series, seasonal.periods = c(24,8760))
plot(tbats_decomp) #decomposição da série
tbats_components <- tbats.components(tbats_decomp)

plot(tbats_decomp$errors) #resíduos
sd_tbats <- sd(tbats_decomp$errors)
tsdisplay(tbats_decomp$errors, lag.max = 50)

Box.test(tbats_decomp$errors, type = "Ljung-Box") #se p-valor <0.05, rejeito H0 (erros iid) com 5% de significância
tbats_dist <- fitdist(as.vector(tbats_decomp$errors),"norm")
plot(tbats_dist)

tbats_signal <- tbats_components[,2] + tbats_components[,3]+tbats_components[,4]

remainder_mc(100,tbats_signal,0,sd_tbats)
MBB_function(100,tbats_signal,tbats_decomp$errors,48,31,T)

# 
# ##
# ##save as png
# ##
# # dev.copy(png, file = "svplots.png", width=10, height=50, units="in", res=400)
# # dev.off()

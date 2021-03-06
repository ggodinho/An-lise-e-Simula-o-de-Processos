#Compara��o entre modelo vigente e simula��es realizadas:

benchmark_comp <- function(series,series_mean,benchmark,benchmark_mean,label){
  plot(series_mean, type = "l",ylim = c(0,1), col="blue",ylab = "Capacity Factor",
       xlab = "Hour",main = paste("Current Model x ",label,sep=""))
  fan(series,fan.col = sequential_hcl,ln=c(5, 95),lw=2)
  lines(benchmark[1,2:25],col="dark red",lw=3)
  lines(benchmark[2,2:25],col="red",lw=3)
  lines(benchmark[3,2:25],col="orange",lw=3)
  axis(side=1, at=c(0,1:24))
  legend("topleft",legend = c("Scenario 1: 29% probability",
                              "Scenario 2: 32% probability",
                              "Scenario 3: 39% probability"),
          col=c("orange","red","dark red"), lwd = c(3,3,3), lty= c(1,1,1), bty="n",
         title = "Benchmark Scenarios", cex = c(0.8,0.8,0.8))
  
  #Comparando Valor Esperado dos Cen�rios
  plot(benchmark_mean, type = "l", ylim = c(0,1), col="red",ylab = "Capacity Factor",
       xlab = "Hour", main = paste("Expected Value: Current Model x ",label,sep=""),
       lwd=2)
  lines(series_mean,col="blue",lw=2,lty=2)
  axis(side=1, at=c(0,1:24))
  legend("bottomright",legend = c("Current Model",label),
         col=c("red","blue"), lwd = c(2,2), lty= c(1,2), bty="n")
  
  # Comparando a densidade de distribui��o das s�ries
  #d1=density(series[1:744])
  d1=density(series_2017[1:744])
  d2=density(series)           
  d3=density(benchmark[,2:25])
  
  plot(d1,ylim=c(0,6), lwd = 1, xlab = "Capacity Factor",xlim=c(0,1),
       main="Density")
  polygon(d1,col=alpha("gray",0.5))
  lines(d2, col = "blue",lty=2,lwd=2)
  lines(d3, col = "red",lwd=2)
  legend("topright", legend = c("Observed Values 2017",label,"Benchmark Model"),
         bty="n", fill=c(alpha("gray",0.5), NA,NA),xpd=T,lty=c(NA,2,1),
         border = c("black", NA,NA),horiz=F,col=c(0,"blue","red"),
         pch = c(0,-1,-1),lwd = c(NA,2,2))
  
}

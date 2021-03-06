
observedvalues_comp <- function(observed_values,model,
                                label = c("Current Model","MSTL+MonteCarlo",
                                          "MSTL+Bootstrap")){
  if (label =="Benchmark Model"){
    plot(observed_values[1,], type = "l",ylim = c(0,1),col=alpha("gray",0.5),xlab="Hour",
         ylab="Capacity Factor",main = paste("Observed Values 2018 x ",label,sep=""))
    for(i in 2:31){
      lines(observed_values[i,],col=alpha("black",0.5))
      }
    lines(model[1,2:25],col="dark red",lw=3)
    lines(model[2,2:25],col="red",lw=3)
    lines(model[3,2:25],col="orange",lw=3)
    axis(side=1, at=c(0,1:24))
  } else {
    plot(observed_values[1,], type = "l",ylim = c(0,1),col=alpha("gray",0.5),xlab="Hour",
         ylab="Capacity Factor",main = paste("Observed Values 2018 x ",label,sep=""),lwd=0.5)
    fan(model, ln=c(5, 25, 50, 75, 95), alpha=0.1,ln.col="yellow",
        lwd= 5,fan.col = colorRampPalette (c("tomato")))
    for(i in 2:31){
      lines(observed_values[i,],col=alpha("black",0.5),lwd=0.5)
    }
    fan(model, ln=c(5, 25, 50, 75, 95), fan.col = colorRampPalette ("tomato"),alpha=0.0,ln.col="yellow",
    lwd= 5)
    axis(side=1, at=c(0,1:24))
  }
}



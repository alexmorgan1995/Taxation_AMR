library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    dS1 = - beta*S1*Is1 - beta*S1*Ir1 + rx*Ir1 + (sigma*rt*Is1 + ((1-sigma)*rw))*Is1   
    dIs1 = beta*S1*Is1 - (sigma*rt*Is1 + ((1-sigma)*rw))*Is1 - (sigma*usr*Is1) + ((1-sigma)*urs*Ir1)  
    dIr1 = beta*S1*Ir1 - rx*Ir1 + (sigma*usr*Is1) - ((1-sigma)*urs*Ir1)  
    
    dS2 = - beta*S2*Is2 - beta*S2*Ir2 + rx*Ir2 + (sigma*rt*Is2 + ((1-sigma)*rw))*Is2   
    dIs2 = beta*S2*Is2 - (sigma*rt*Is2 + ((1-sigma)*rw))*Is2 - (sigma*usr*Is2) + ((1-sigma)*urs*Ir2)  
    dIr2 = beta*S2*Ir2 - rx*Ir2 + (sigma*usr*Is2) - ((1-sigma)*urs*Ir2)  
    
    dS3 = - beta*S3*Is3 - beta*S3*Ir3 + rx*Ir3 + (sigma*rt*Is3 + ((1-sigma)*rw))*Is3   
    dIs3 = beta*S3*Is3 - (sigma*rt*Is3 + ((1-sigma)*rw))*Is3 - (sigma*usr*Is3) + ((1-sigma)*urs*Ir3)  
    dIr3 = beta*S3*Ir3 - rx*Ir3 + (sigma*usr*Is3) - ((1-sigma)*urs*Ir3) 
    
    CumIs1 = beta*S1*Is1
    CumIr1 = beta*S1*Ir1
    
    CumIs2 = beta*S2*Is2
    CumIr2 = beta*S2*Ir2
    
    CumIs3 = beta*S3*Is3
    CumIr3 = beta*S3*Ir3
    
    return(list(c(dS1,dIs1,dIr1,
                  dS2,dIs2,dIr2,
                  dS3,dIs3,dIr3), 
                CumIs1, CumIr1,
                CumIs2, CumIr2,
                CumIs3, CumIr3))
  })
}

# Parameters and Initial Conditions ---------------------------------------

init <- c(S1 = 0.99, Is1 = (1-0.99), Ir1 = 0,
          S2 = 0.99, Is2 = (1-0.99), Ir2 = 0,
          S3 = 0.99, Is3 = (1-0.99), Ir3 = 0)

parms = c(beta = 0.20, 
          sigma = 0.5,
          rw = 1/14,
          rx = 1/12,
          rt = 1/7,
          usr = 0.01,
          urs = 0.01)

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms))

ggplot(data = testrun, aes(x = time, y = Ir1/(Ir1+Is1))) + geom_line()

ggplot(data = testrun, aes(x = time, y = Ir2/(Ir2+Is2))) + geom_line()

ggplot(data = testrun, aes(x = time, y = Ir3/(Ir3+Is3))) + geom_line()

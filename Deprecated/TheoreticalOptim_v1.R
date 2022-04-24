library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    if(t > t_n) {
    sigma <- sigma*0.2  
    } 
    
    dS = - beta*S*X - beta*S*Ir1 - beta*S*Ir2 + rx1*Ir1 + rx2*Ir2 + (sigma*rt + (1-sigma)*rw)*X
    dIr1 = beta*S*Ir1 - rx1*Ir1 - m_rx*Ir1 + m_xr*X
    dIr2 = beta*S*Ir2 - rx2*Ir2 - m_rx*Ir2 + m_xr*X
    dX = beta*S*X - (sigma*rt + (1-sigma)*rw)*X + m_rx*Ir1 - m_xr*X + m_rx*Ir2 - m_xr*X
  
    return(list(c(dS,dIr1,dIr2,dX)))
  })
}

# Parameters and Initial Conditions ---------------------------------------

init <- c(S = 0.99, Ir1 = 0, Ir2 =0, X = 1-0.99)

parms = c(beta = 0.20, 
          sigma1 = 0.5,
          sigma2 = 0.5,
          rw = 1/14,
          rx1 = 1/12,
          rx2 = 1/12,
          rt = 1/7,
          m_xr = 0.01,
          m_rx = 0.01,
          t_n = 100)

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms))

testrun[6:8] <- testrun[,3:5]/rowSums(testrun[,3:5])

comb_prop = melt(testrun, id.vars = "time", measure.vars = tail(colnames(testrun), 3))

ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line()


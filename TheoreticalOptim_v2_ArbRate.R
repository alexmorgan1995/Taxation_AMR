library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    if(t > t_n) {
    sigma_1 <- sigma_1*0.2  
    sigma_2 <- sigma_2*0.2  
    }
    
    dS = - beta*S*(X+R1+R2) + (1-sigma_1+sigma_2)*mu_x*X + (sigma_1+sigma_2)*mu_t*X +
      (1-sigma_2)*mu_r1*R1 + sigma_2*mu_t*R1 +
      (1-sigma_1)*mu_r2*R2 + sigma_1*mu_t*R2
    
    dX = beta*S*X - (1-sigma_1+sigma_2)*mu_x*X - (sigma_1+sigma_2)*mu_t*X + eta_rx_1*R1 + eta_rx_2*R2 - 
      eta_xr_1*X - eta_xr_2*X
    
    dR1 = beta*S*R1 - (1-sigma_2)*mu_r1*R1 - sigma_2*mu_t*R1 - eta_rx_1*R1 + eta_xr_1*X
    
    dR2 = beta*S*R2 - (1-sigma_1)*mu_r2*R2 - sigma_1*mu_t*R2 - eta_rx_2*R2 + eta_xr_2*X
    
    return(list(c(dS,dX,dR1,dR2)))
  })
}

# Parameters and Initial Conditions ---------------------------------------

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0)

parms = c(beta = 0.1, 
          sigma_1 = 0.4,
          sigma_2 = 0.4,
          mu_x = 1/20,
          mu_r1 = 1/10,
          mu_r2 = 1/10,
          mu_t = 1/4,
          eta_xr_1 = 0.01,
          eta_rx_1 = 0.01,
          eta_xr_2 = 0.05,
          eta_rx_2 = 0.01,
          t_n = 100)

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms))
testrun[6:8] <- testrun[,3:5]/rowSums(testrun[,3:5])
testrun$avgRes <- rowMeans(testrun[,7:8])
colnames(testrun)[6:9] <- c("WT", "Res1", "Res2", "AvgRes")
comb_prop = melt(testrun, id.vars = "time", measure.vars = tail(colnames(testrun), 4))

ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line() + theme_bw() +
  scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0)) +
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") 

library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    if(t > t_n) {
      sigma_1 <- sigma_1*eff_tax1  
      sigma_2 <- sigma_2*eff_tax2
    }
    
    dS = - beta*S*(X+R1+R2) + (1-sigma_1+sigma_2)*mu_x*X + (sigma_1+sigma_2)*mu_t*X +
      (1-sigma_2)*mu_r1*R1 + sigma_2*mu_t*R1 +
      (1-sigma_1)*mu_r2*R2 + sigma_1*mu_t*R2
    
    dX = beta*S*X - (1-sigma_1+sigma_2)*mu_x*X - (sigma_1+sigma_2)*mu_t*X + 
      (1-sigma_1+sigma_2)*eta_rx_1*R1 + (1-sigma_1+sigma_2)*eta_rx_2*R2 - 
      sigma_1*eta_xr_1*X - sigma_2*eta_xr_2*X
    
    dR1 = beta*S*R1 - (1-sigma_2)*mu_r1*R1 - sigma_2*mu_t*R1 - 
      (1-sigma_1+sigma_2)*eta_rx_1*R1 + sigma_1*eta_xr_1*X
    
    dR2 = beta*S*R2 - (1-sigma_1)*mu_r2*R2 - sigma_1*mu_t*R2 - 
      (1-sigma_1+sigma_2)*eta_rx_2*R2 + sigma_2*eta_xr_2*X
    
    dC_X = 0
    dC_R1 = 0
    dC_R2 = 0
    
    #Calculating the Proportion Integrals
    if(t > t_n) {
      dC_X = (beta*S*X + eta_rx_1*R1 + eta_rx_2*R2)/((beta*S*X + eta_rx_1*R1 + eta_rx_2*R2)+
                                                       beta*S*R1 + eta_xr_1*X + beta*S*R2 + eta_xr_2*X)
      dC_R1 = beta*S*R1 + eta_xr_1*X/((beta*S*X + eta_rx_1*R1 + eta_rx_2*R2)+
                                        beta*S*R1 + eta_xr_1*X + beta*S*R2 + eta_xr_2*X)
      dC_R2 = beta*S*R2 + eta_xr_2*X/((beta*S*X + eta_rx_1*R1 + eta_rx_2*R2)+
                                        beta*S*R1 + eta_xr_1*X + beta*S*R2 + eta_xr_2*X)
    }
    
    return(list(c(dS,dX,dR1,dR2,dC_X,dC_R1,dC_R2)))
  })
}

# Parameters and Initial Conditions ---------------------------------------

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0,
          C_X = 0, C_R1 = 0, C_R2 =0)

parms = c(beta = 0.3, 
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
          eff_tax1 = 1,
          eff_tax2 = 0.2,
          t_n = 100)

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms))

ggplot(testrun, aes(time, R1)) + geom_line()

testrun[9:11] <- testrun[,3:5]/rowSums(testrun[,3:5])
testrun$avgRes <- rowMeans(testrun[,10:11])
colnames(testrun)[9:12] <- c("WT", "Res1", "Res2", "AvgRes")
comb_prop = melt(testrun, id.vars = "time", measure.vars = tail(colnames(testrun), 4))

ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line() + theme_bw() +
  scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0)) +
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") 
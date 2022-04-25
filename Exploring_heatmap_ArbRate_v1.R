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
    
    dX = beta*S*X - (1-sigma_1+sigma_2)*mu_x*X - (sigma_1+sigma_2)*mu_t*X + eta_rx_1*R1 + eta_rx_2*R2 - 
      eta_xr_1*X - eta_xr_2*X
    
    dR1 = beta*S*R1 - (1-sigma_2)*mu_r1*R1 - sigma_2*mu_t*R1 - eta_rx_1*R1 + eta_xr_1*X
    
    dR2 = beta*S*R2 - (1-sigma_1)*mu_r2*R2 - sigma_1*mu_t*R2 - eta_rx_2*R2 + eta_xr_2*X
    
    #Calculating the Proportion Integrals
    
    return(list(c(dS,dX,dR1,dR2)))
  })
}

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:5]),
                      sum(rowMeans(data_temp[4:5])),
                      sum(data_temp[4]),
                      sum(data_temp[5])),5)
  return(out_vec)
}

# 3x3 Trajectory Plots ----------------------------------------------------

comb_tax <- expand.grid(efftax1 = seq(0.25,0.75,0.25),
                        efftax2 = seq(0.25,0.75,0.25))

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0)

parms1 = c(beta = 0.5, sigma_1 = 0.35, sigma_2 = 0.35, mu_x = 1/15, mu_r1 = 1/10,
          mu_r2 = 1/12, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
          eta_xr_2 = 0.01, eta_rx_2 = 0.01, eff_tax1 = .2,
          eff_tax2 = .2, t_n = 500)

data_store <- list()

for(i in 1:nrow(comb_tax)){
  parms1[c("eff_tax1", "eff_tax2")] <- comb_tax[i,]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 1000), parms = parms1))
  
  run$AvgRes <- rowMeans(run[,4:5])
  run$TotInf <- rowSums(run[,3:5])
  
  comb_prop = melt(run, id.vars = "time", measure.vars = c("X","R1","R2", "AvgRes", "TotInf"))
  
  data_store[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
                      scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.9)) +
                      theme(legend.text=element_text(size=12), legend.title = element_text(size=12), axis.text=element_text(size=12), plot.title = element_text(size=11),
                      axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
                      legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + scale_linetype_manual(values = c("solid", "solid","solid", "longdash", "longdash"), guide = 'none') +
    ggtitle(paste0("Tax 1: ", (1-parms1[["eff_tax1"]])*100, "% | Tax2: ", (1-parms1[["eff_tax2"]])*100, "%")) 
}

ggarrange(plotlist = data_store, ncol = 3, nrow = 3, 
          common.legend = T, legend = "bottom")

# Mutation Rate Exploration -----------------------------------------------

explore_vec <- seq(0, 1, 0.25)

for(i in 1:length(explore_vec)) {
  
  
}



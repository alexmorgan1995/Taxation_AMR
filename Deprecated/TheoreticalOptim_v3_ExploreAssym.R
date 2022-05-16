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
    
    dC_X = 0
    dC_R1 = 0
    dC_R2 = 0
    
    #Calculating the Proportion Integrals
    if(t > t_n) {
      dC_X = beta*S*X + eta_rx_1*R1 + eta_rx_2*R2
      dC_R1 = beta*S*R1 + eta_xr_1*X
      dC_R2 = beta*S*R2 + eta_xr_2*X
    }
    
    return(list(c(dS,dX,dR1,dR2,dC_X,dC_R1,dC_R2)))
  })
}

# Parameters and Initial Conditions ---------------------------------------

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0,
          C_X = 0, C_R1 = 0, C_R2 =0)

parms = c(beta = 0.35, 
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
          eff_tax2 = 1,
          t_n = 100)

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms))

ggplot(testrun, aes(time, R1)) + geom_line()

testrun$avgRes <- rowMeans(testrun[,4:5]); colnames(testrun)[9] <- c("AvgRes")

comb_prop = melt(testrun, id.vars = "time", measure.vars = c("X","R1","R2","AvgRes"))

ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line() + theme_bw() +
  scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.7)) +
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom")

# Exploring Assymetry -----------------------------------------------------
# Ranging Beta ------------------------------------------------------------
#Altering Beta

beta_range <- c(0.1, 0.2, 0.3, 0.4)
parms1 = c(beta = 0.35, sigma_1 = 0.4, sigma_2 = 0.4, mu_x = 1/20, mu_r1 = 1/10,
           mu_r2 = 1/10, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
           eta_xr_2 = 0.01, eta_rx_2 = 0.05, eff_tax1 = 0.2,
           eff_tax2 = 0.2, t_n = 100)

data_store <- list()

for(i in 1:length(beta_range)){
  parms1["beta"] <- beta_range[i]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms1))
  run$AvgRes <- rowMeans(run[,4:5]); run$TotInf <- rowSums(run[,3:5])
  
  comb_prop = melt(run, id.vars = "time", measure.vars = c("X","R1","R2", "AvgRes", "TotInf"))
  
  data_store[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.8)) +
    theme(legend.text=element_text(size=12), legend.title = element_text(size=12), axis.text=element_text(size=12), plot.title = element_text(size=12),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + scale_linetype_manual(values = c("solid", "solid","solid", "longdash", "longdash"), guide = 'none') +
    ggtitle(paste0("Beta: ", parms1["beta"])) 
}

#Plotting Scenarios

ggarrange(plotlist = data_store, ncol = 4, nrow = 1, common.legend = T, legend = "bottom")

# Ranging Eta ------------------------------------------------------------

parms1 = c(beta = 0.35, sigma_1 = 0.4, sigma_2 = 0.4, mu_x = 1/20, mu_r1 = 1/10,
           mu_r2 = 1/10, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
           eta_xr_2 = 0.05, eta_rx_2 = 0.01, eff_tax1 = 0.2,
           eff_tax2 = 0.2, t_n = 100)

eta_diff <- parms1 

data_store_eta <- list()

for(i in 1:3){
  eta_diff[c("eta_xr_2", "eta_rx_2")] <- list(c(0.01, 0.01),
                                              c(0.05, 0.01),
                                              c(0.1, 0.01))[[i]]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = eta_diff))
  run$AvgRes <- rowMeans(run[,4:5]); run$TotInf <- rowSums(run[,3:5])
  
  comb_prop = melt(run, id.vars = "time", measure.vars = c("X","R1","R2", "AvgRes", "TotInf"))
  
  data_store_eta[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.8)) +
    theme(legend.text=element_text(size=12), legend.title = element_text(size=12), axis.text=element_text(size=12), plot.title = element_text(size=12),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + scale_linetype_manual(values = c("solid", "solid","solid", "longdash", "longdash"), guide = 'none') +
    ggtitle(paste0("eta_xr2: ", eta_diff["eta_xr_2"], " | eta_rx2: ", eta_diff["eta_rx_2"])) 
}

ggarrange(plotlist = data_store_eta, ncol = 3, nrow = 1, common.legend = T, legend = "bottom")

# Differences in Recovery Rate --------------------------------------------

parms1 = c(beta = 0.35, sigma_1 = 0.4, sigma_2 = 0.4, mu_x = 1/20, mu_r1 = 1/10,
           mu_r2 = 1/10, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
           eta_xr_2 = 0.01, eta_rx_2 = 0.01, eff_tax1 = 0.2,
           eff_tax2 = 0.2, t_n = 100)

mu_diff <- parms1 

data_store_mu <- list()

for(i in 1:2){
  mu_diff[c("mu_r1", "mu_r2")] <- list(c(1/10, 1/10),
                                              c(1/8, 1/15))[[i]]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = mu_diff))
  run$AvgRes <- rowMeans(run[,4:5]); run$TotInf <- rowSums(run[,3:5])
  
  comb_prop = melt(run, id.vars = "time", measure.vars = c("X","R1","R2", "AvgRes", "TotInf"))
  
  data_store_mu[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.8)) +
    theme(legend.text=element_text(size=12), legend.title = element_text(size=12), axis.text=element_text(size=12), plot.title = element_text(size=12),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + scale_linetype_manual(values = c("solid", "solid","solid", "longdash", "longdash"), guide = 'none') +
    ggtitle(paste0("mu_r1: ", mu_diff["mu_r1"], " | mu_r2: ", mu_diff["mu_r2"])) 
}

ggarrange(plotlist = data_store_mu, ncol = 2, nrow = 1, common.legend = T, legend = "bottom")


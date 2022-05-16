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
    
    dX = - beta*X*(W+(R1*c1)+(R2*c2)) + (1-sigma_1+sigma_2)*mu_w*W + (sigma_1+sigma_2)*mu_t*W +
      (1-sigma_2)*mu_r*R1 + sigma_2*mu_t*R1 +
      (1-sigma_1)*mu_r*R2 + sigma_1*mu_t*R2
    
    dW = beta*X*W - (1-sigma_1+sigma_2)*mu_w*W - (sigma_1+sigma_2)*mu_t*W + 
      (1-sigma_1+sigma_2)*eta_rw*c1*R1 + (1-sigma_1+sigma_2)*eta_rw*c2*R2 - 
      sigma_1*eta_wr*W - sigma_2*eta_wr*W
    
    dR1 = beta*X*R1*c1 - (1-sigma_2)*mu_r*R1 - sigma_2*mu_t*R1 - (1-sigma_1+sigma_2)*eta_rw*c1*R1 + sigma_1*eta_wr*W
    
    dR2 = beta*X*R2*c2 - (1-sigma_1)*mu_r*R2 - sigma_1*mu_t*R2 - (1-sigma_1+sigma_2)*eta_rw*c2*R2 + sigma_2*eta_wr*W
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX,dW,dR1,dR2)))
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

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0)

parms1 = c(beta = 1, sigma_1 = 0.35, sigma_2 = 0.35, mu_w = 1/15, mu_r = 1/10, 
           mu_t = 1/5, eta_wr = 0.02, eta_rw= 0.02, c1= 0.9, c2 = 0.85,
           eff_tax1 = 1, eff_tax2 = 0.5, t_n = 1000)

data_store <- list()

for(i in 1:nrow(comb_tax)){
  parms1[c("eff_tax1", "eff_tax2")] <- comb_tax[i,]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms1))
  
  run$AvgRes <- rowMeans(run[,4:5])
  run$TotInf <- rowSums(run[,3:5])
  
  comb_prop = melt(run, id.vars = "time", measure.vars = c("X","R1","R2", "AvgRes", "TotInf"))
  
  data_store[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
                      scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,1)) +
                      theme(legend.text=element_text(size=12), legend.title = element_text(size=12), axis.text=element_text(size=12), plot.title = element_text(size=11),
                      axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
                      legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + scale_linetype_manual(values = c("solid", "solid","solid", "dashed", "dashed"), guide = 'none') +
    ggtitle(paste0("Tax 1: ", (1-parms1[["eff_tax1"]])*100, "% | Tax2: ", (1-parms1[["eff_tax2"]])*100, "%")) 
}

ggarrange(plotlist = data_store, ncol = 3, nrow = 3, 
          common.legend = T, legend = "bottom")

# Mutation Rate Exploration -----------------------------------------------

explore_vec <- seq(0, 1, 0.25)

for(i in 1:length(explore_vec)) {
  
  
}



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

# Parameters and Initial Conditions ---------------------------------------

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0)

parms1 = c(beta = 1, sigma_1 = 0.35, sigma_2 = 0.35, mu_w = 1/15, mu_r = 1/10, 
           mu_t = 1/5, eta_wr = 0.02, eta_rw= 0.02, c1= 0.9, c2 = 0.85,
           eff_tax1 = 0.5, eff_tax2 = 1, t_n = 1000)

# Run the Model - R1 Reductions -----------------------------------------------------------

R1_usage <- seq(0.6, 1, by = 0.1)

p_list <- list()

for(i in 1:length(R1_usage)) {
  parms1["eff_tax1"] <- R1_usage[i]

  run <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms1))
  
  run$AvgRes <- rowMeans(run[,4:5])
  run$TotInf <- rowSums(run[,3:5])
  colnames(run)[2:7] <- c("X", "W", "High Res", "Low Res", "Average", "TotInf")
  comb_prop = melt(run, id.vars = "time", measure.vars =  c("W", "High Res", "Low Res", "TotInf"))
  
  p_list[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,1)) +
    theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), plot.title = element_text(size=11),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + scale_linetype_manual(values = c("solid","solid","solid", "dashed"), guide = 'none') +
    ggtitle(paste0("Reductions to HighRes: ", (1-round(parms1[["eff_tax1"]], 2))*100, "%")) + 
    annotate("rect", xmin=1000, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="red")
}
  
ggarrange(plotlist = p_list, nrow =1, ncol = length(R1_usage), common.legend = T, legend = "bottom")

# Run the Model - R2 Reductions -----------------------------------------------------------

R2_usage <- seq(0.6, 1, by = 0.1)
p_list <- list()

for(i in 1:length(R1_usage)) {
  parms1["eff_tax2"] <- R2_usage[i]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms1))
  
  run$AvgRes <- rowMeans(run[,4:5])
  run$TotInf <- rowSums(run[,3:5])
  colnames(run)[2:7] <- c("X", "W", "High Res", "Low Res", "Average", "TotInf")
  comb_prop = melt(run, id.vars = "time", measure.vars =  c("W", "High Res", "Low Res", "TotInf"))
  
  p_list[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,1)) +
    theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), plot.title = element_text(size=11),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + scale_linetype_manual(values = c("solid","solid","solid", "dashed"), guide = 'none') +
    ggtitle(paste0("Reductions to LowRes: ", (1-round(parms1[["eff_tax2"]], 2))*100, "%")) + 
    annotate("rect", xmin=1000, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="red")
}

ggarrange(plotlist = p_list, nrow =1, ncol = length(R1_usage), common.legend = T, legend = "bottom")
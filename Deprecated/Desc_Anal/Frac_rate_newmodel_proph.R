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
    
    dPr1 = 
    
    dPr2 = 
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX,dW,dR1,dR2)))
  })
}

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:5]),
                      sum(data_temp[4:5]),
                      sum(data_temp[3]),
                      sum(rowMeans(data_temp[4:5])),
                      sum(data_temp[4]),
                      sum(data_temp[5])),5)
  return(out_vec)
}

# Parameters and Initial Conditions ---------------------------------------

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0)

parms = c(beta = 1, sigma_1 = 0.35, sigma_2 = 0.35, mu_w = 1/15, mu_r = 1/10, 
          mu_t = 1/5, eta_wr = 0.02, eta_rw= 0.02, c1= 0.9, c2 = 0.85,
          eff_tax1 = 1, eff_tax2 = 0.5, t_n = 1000)

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms))

testrun$avgRes <- rowMeans(testrun[,4:5]); colnames(testrun)[6] <- c("AvgRes")
colnames(testrun)[2:6] <- c("X", "W", "High Res", "Low Res", "Average")
comb_prop = melt(testrun, id.vars = "time", measure.vars = colnames(testrun)[2:6])

ggplot(data = comb_prop, aes(x = time, y = value, color = variable, size = variable)) + geom_line() + theme_bw() +
  annotate("rect", xmin=1000, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="red") + 
  scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,1)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + 
  scale_size_manual(values = c(1,1,2,2,2))


# Plotting the Trajectory Plots -------------------------------------------

tax_discrete <- list(c(1, 1),
                     c(0.5, 0.5),
                     c(0.5, 1),
                     c(1, 0.5))

parms1 = c(beta = 1, sigma_1 = 0.35, sigma_2 = 0.35, mu_w = 1/15, mu_r = 1/10, 
          mu_t = 1/5, eta_wr = 0.02, eta_rw= 0.02, c1= 0.9, c2 = 0.85,
          eff_tax1 = 1, eff_tax2 = 0.5, t_n = 1000)

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0)

p_data <- list()

for(i in 1:4) {
  parms1[c("eff_tax1","eff_tax2")] <- tax_discrete[[i]]
  testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms1))
  
  testrun$AvgRes <- rowMeans(testrun[,4:5])
  colnames(testrun)[2:6] <- c("X", "W", "High Res", "Low Res", "Average")

  comb_prop = melt(testrun, id.vars = "time", measure.vars = colnames(testrun)[2:6])
  
  p_data[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable, size = variable)) + geom_line() + theme_bw() + 
    annotate("rect", xmin=1000, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="red") +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.9)) +
    theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), plot.title = element_text(size=11),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + 
    ggtitle(c("No Tax",
              "Both Taxed",
              "Only High Res Taxed",
              "Only Low Res Taxed")[i]) + 
    scale_size_manual(values = c(1,1,2,2,2)) 
  
  print(integral(testrun, parms1["t_n"]))
}

ggarrange(plotlist = p_data, common.legend = T, ncol = 4, nrow = 1, legend = "bottom")

# Bivariate Optimisation - Tax --------------------------------------------------

#Reductions to Fraction Treated 1  (Sigma 1) and Fraction Treated 2 (Sigma 2)

comb_tax <- expand.grid(efftax1 = seq(0,1,0.1),
                          efftax2 = seq(0,1,0.1))

parms1 = c(beta = 1, sigma_1 = 0.35, sigma_2 = 0.35, mu_w = 1/15, mu_r = 1/10, 
           mu_t = 1/5, eta_wr = 0.02, eta_rw= 0.02, c1= 0.9, c2 = 0.85,
           eff_tax1 = 1, eff_tax2 = 0.5, t_n = 1000)

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0)

data_store <- data.frame(matrix(ncol = 6, nrow = length(comb_tax)))

for(i in 1:nrow(comb_tax)){
  parms1[c("eff_tax1", "eff_tax2")] <- comb_tax[i,]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms1))
  data_store[i,] <- c((1-comb_tax[i,])*100,
                      integral(run, parms1["t_n"]))
}

colnames(data_store) <- c("class1", "class2", "wt_res", "avgres", "res1", "res2")

#Plotting Scenarios

plot_list <- list()

for(i in 1:4){
  plot_list[[i]] <- ggplot(data_store, aes_string(x = "class1", y = "class2", 
                                           fill = c("wt_res", "avgres", "res1", "res2")[i])) + 
    geom_tile() + labs(fill = c("Tot_Inf", "AvgRes", "Res1", "Res 2")[i]) +
    scale_x_continuous(name = bquote("Usage Reduction High Res Class - R1 (%)"), expand = c(0, 0)) +  
    scale_y_continuous(name = bquote("Usage Reduction Low Res Class - R2 (%)"), expand = c(0, 0)) +
    scale_fill_viridis(direction = -1) 
}

ggarrange(plotlist = plot_list)

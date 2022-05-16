library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    if(t > t_n | switch == 1) {
    #h_res <- eval(as.name(c("R1", "R2", "R3")[which.max(c(R1, R2, R3))]))
    #m_res <- eval(as.name(c("R1", "R2", "R3")[setdiff(1:3, c(which.min(c(R1, R2, R3)), which.max(c(R1, R2, R3))))]))
    #l_res <- eval(as.name(c("R1", "R2", "R3")[which.min(c(R1, R2, R3))]))
      
    sigma_1 <- sigma_1*eff_tax1 
    sigma_2 <- sigma_2*eff_tax2
    sigma_3 <- sigma_3*eff_tax3
    }
    
    if(t > t_n | switch == 2) {
      sigma_1 <- sigma_1*eff_tax 
      sigma_2 <- sigma_2*eff_tax
      sigma_3 <- sigma_2*eff_tax3
    }
    
    
    dX = - beta*X*(W+(R1)+(R2)+(R3)) + (1-sigma_1+sigma_2+sigma_3)*mu_w*W + (sigma_1+sigma_2+sigma_3)*mu_t*W +
      (1-sigma_2+sigma_3)*mu_r1*R1 + (sigma_2+sigma_3)*mu_t*R1 +
      (1-sigma_1+sigma_3)*mu_r2*R2 + (sigma_1+sigma_3)*mu_t*R2 + 
      (1-sigma_1+sigma_2)*mu_r3*R3 + (sigma_1+sigma_2)*mu_t*R3
    
    dW = beta*X*W - (1-sigma_1+sigma_2+sigma_3)*mu_w*W - (sigma_1+sigma_2+sigma_3)*mu_t*W + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R1 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R2 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R3 - 
      sigma_1*eta_wr*W - sigma_2*eta_wr*W - sigma_3*eta_wr*W
    
    dR1 = beta*X*R1 - (1-sigma_2+sigma_3)*mu_r1*R1 - (sigma_2+sigma_3)*mu_t*R1 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R1 + sigma_1*eta_wr*W
    
    dR2 = beta*X*R2 - (1-sigma_1+sigma_3)*mu_r2*R2 - (sigma_1+sigma_3)*mu_t*R2 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R2 + sigma_2*eta_wr*W
    
    dR3 = beta*X*R3 - (1-sigma_1+sigma_2)*mu_r3*R3 - (sigma_1+sigma_2)*mu_t*R3 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R3 + sigma_3*eta_wr*W
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX,dW,dR1,dR2,dR3)))
  })
}


# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:6]),
                      sum(data_temp[4:6]),
                      sum(data_temp[3]),
                      sum(rowMeans(data_temp[4:6])),
                      sum(data_temp[4]),
                      sum(data_temp[5]),
                      sum(data_temp[6])),5)
  return(out_vec)
}

# Parameters and Initial Conditions ---------------------------------------

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 = 0, R3 = 0)

parms = c(beta = 5, sigma_1 = 0.25, sigma_2 = 0.25, sigma_3 = 0.25,
          mu_w = 1/14, mu_r1 = 1/13, mu_r2 = 1/12.75, mu_r3 = 1/12, 
          mu_t = 1/7, eta_wr = 0.01, eta_rw= 0.01, 
          eff_tax1 = 0.5, eff_tax2 = 1, eff_tax3 = 1, 
          t_n = 3000)

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 5000), parms = parms))

testrun$AvgRes <- rowMeans(testrun[,4:6]);

sum(tail(testrun, 1)[2:6])

comb_prop = melt(testrun, id.vars = "time", measure.vars = c("X","W","R1","R2","R3","AvgRes"))

integral(testrun,3000)

ggplot(data = comb_prop, aes(x = time, y = value, color = variable, lty = variable)) + geom_line() + theme_bw() +
  annotate("rect", xmin=3000, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="red") + 
  scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,1)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") +
  scale_linetype_manual(values = c("W" = "dashed","R1" = "solid","R2" = "solid",
                                   "R3" = "solid","AvgRes" = "dashed"), guide = "none") + 
  guides(colour = guide_legend(show = FALSE))

# Plotting the Trajectory Plots -------------------------------------------

tax_discrete <- as.matrix(expand.grid(c(1,0.5),
                                      c(1,0.5),
                                      c(1,0.5)))

parms1 = c(beta = 5, sigma_1 = 0.25, sigma_2 = 0.25,sigma_3 = 0.25,
          mu_w = 1/12, mu_r = 1/10, 
          mu_t = 1/7, eta_wr = 0.01, eta_rw= 0.01, 
          c1= 0.95, c2 = 0.94, c3 = 0.91,
          eff_tax1 = 0.5, eff_tax2 = 0.5, eff_tax3 = 0.5, 
          t_n = 3000)

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0, R3 =0)

p_data <- list()

for(i in 1:nrow(tax_discrete)) {
  parms1[c("eff_tax1","eff_tax2","eff_tax3")] <- tax_discrete[i,]
  
  testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 7500), parms = parms1))
  
  testrun$AvgRes <- rowMeans(testrun[,4:6])
  comb_prop = melt(testrun, id.vars = "time", measure.vars = c("W","R1","R2","R3","AvgRes"))
  
  p_data[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line() + theme_bw() +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.9)) +
    theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), plot.title = element_text(size=11),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + 
    ggtitle(paste0("Tax1: ", (1-parms1["eff_tax1"])*100, "%", " | Tax2: ", (1-parms1["eff_tax2"])*100, "%",
                   " | Tax3: ", (1-parms1["eff_tax3"])*100, "%")) 
  
  print(paste0("Tax1: ", (1-parms1["eff_tax1"])*100, "%", " | Tax2: ", (1-parms1["eff_tax2"])*100, "%",
               " | Tax3: ", (1-parms1["eff_tax3"])*100, "%"))
  print(integral(testrun, parms1["t_n"]))
}

ggarrange(plotlist = p_data, common.legend = T, ncol = 4, nrow = 2, legend = "bottom")

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
  data_store[i,] <- c(1-comb_tax[i,],
                      integral(run, parms1["t_n"]))
}

colnames(data_store) <- c("class1", "class2", "wt_res", "avgres", "res1", "res2")

#Plotting Scenarios

plot_list <- list()

for(i in 1:4){
  plot_list[[i]] <- ggplot(data_store, aes_string(x = "class1", y = "class2", 
                                           fill = c("wt_res", "avgres", "res1", "res2")[i])) + 
    geom_tile() + labs(fill = c("Tot_Inf", "AvgRes", "Res1", "Res 2")[i]) +
    scale_x_continuous(name = bquote("Usage Reduction Class 1 (%)"), expand = c(0, 0)) +  
    scale_y_continuous(name = bquote("Usage Reduction Class 2 (%)"), expand = c(0, 0)) +
    scale_fill_viridis(direction = -1) 
}

ggarrange(plotlist = plot_list)

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

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms))

testrun$avgRes <- rowMeans(testrun[,4:5]); colnames(testrun)[9] <- c("AvgRes")

comb_prop = melt(testrun, id.vars = "time", measure.vars = c("X","R1","R2","AvgRes"))

ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line() + theme_bw() +
  scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.7)) +
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom")

# Calculating Integrals ---------------------------------------------------

# The Integral of the proportion infecteds are simply the cumulative number of entries into the compartment
# over time

#These last two are what we are trying to integrate over time

Int_cum <- data.frame("time" = testrun$time,
                      "cum_resprop" = rowMeans(testrun[,7:8])/2,
                      "cum_sum" = rowSums(testrun[,6:8]))

# Plotting the Trajectory Plots -------------------------------------------

tax_discrete <- list(c(1, 1),
                     c(0.2, 0.2),
                     c(0.2, 1),
                     c(1, 0.2))

parms1 = c(beta = 0.35, sigma_1 = 0.4, sigma_2 = 0.4, mu_x = 1/20, mu_r1 = 1/10,
           mu_r2 = 1/12, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
           eta_xr_2 = 0.01, eta_rx_2 = 0.01, eff_tax1 = 1,
           eff_tax2 = 1, t_n = 500)

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0,
          C_X = 0, C_R1 = 0, C_R2 =0)

p_data <- list()

for(i in 1:4) {
  parms1[c("eff_tax1","eff_tax2")] <- tax_discrete[[i]]
  testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 1000), parms = parms1))
  
  testrun$avgRes <- rowMeans(testrun[,4:5]); colnames(testrun)[9] <- c("AvgRes")
  comb_prop = melt(testrun, id.vars = "time", measure.vars = c("X","R1","R2","AvgRes"))
  
  p_data[[i]] <- ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line() + theme_bw() +
    scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,0.75)) +
    theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), plot.title = element_text(size=11),
          axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + 
    ggtitle(paste0("Tax1: ", (1-parms1["eff_tax1"])*100, "%", " | Tax2: ", (1-parms1["eff_tax2"])*100, "%")) 
  
  print(paste0(signif(c(sum(tail(testrun[,6:8],1)),
                        rowMeans(tail(testrun[,7:8],1)),
                        tail(testrun[,7],1),
                        tail(testrun[,8],1)))), 5)
}

ggarrange(plotlist = p_data, common.legend = T, ncol = 4, nrow = 1, legend = "bottom")


# Bivariate Optimisation - Tax --------------------------------------------------

#Reductions to Fraction Treated 1  (Sigma 1) and Fraction Treated 2 (Sigma 2)

comb_tax <- expand.grid(efftax1 = seq(0,1,0.1),
                          efftax2 = seq(0,1,0.1))

parms1 = c(beta = 1, sigma_1 = 0.4, sigma_2 = 0.4, mu_x = 1/20, mu_r1 = 1/10,
          mu_r2 = 1/12, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
          eta_xr_2 = 0.01, eta_rx_2 = 0.01, eff_tax1 = 0.2,
          eff_tax2 = 0.2, t_n = 500)

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0,
          C_X = 0, C_R1 = 0, C_R2 =0)

data_store <- data.frame(matrix(ncol = 6, nrow = length(comb_tax)))

for(i in 1:nrow(comb_tax)){
  parms1[c("eff_tax1", "eff_tax2")] <- comb_tax[i,]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 1000), parms = parms1))
  
  run$avgRes <- rowMeans(run[,4:5])
  colnames(run)[9] <- c("AvgPropRes")
  run$resprop <- rowMeans(run[,7:8])
  run$cumsum <- rowSums(run[,6:8])
  
  data_store[i,] <- c((1-comb_tax[i,])*100,
                   tail(run$cumsum,1),
                   tail(run$resprop,1),
                   tail(run[,7],1),
                   tail(run[,8],1)) 
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

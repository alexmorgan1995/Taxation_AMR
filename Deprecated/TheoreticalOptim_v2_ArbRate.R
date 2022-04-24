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
          eff_tax1 = 0.2,
          eff_tax2 = 1,
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

# Calculating Integrals ---------------------------------------------------

# The Integral of the proportion infecteds are simply the cumulative number of entries into the compartment
# over time

#These last two are what we are trying to integrate over time

Int_cum <- data.frame("time" = testrun$time,
                      "cum_resprop" = rowMeans(testrun[,7:8])/2,
                      "cum_sum" = rowSums(testrun[,6:8]))

# Bivariate Optimisation --------------------------------------------------

#Reductions to Fraction Treated 1  (Sigma 1) and Fraction Treated 2 (Sigma 2)

comb_sigma <- expand.grid(sigma_1 = seq(0,1,0.1),
                          sigma_2 = seq(0,1,0.1))

parms1 = c(beta = 0.1, sigma_1 = 0.4, sigma_2 = 0.4, mu_x = 1/20, mu_r1 = 1/10,
          mu_r2 = 1/10, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
          eta_xr_2 = 0.05, eta_rx_2 = 0.01, eff_tax1 = 0.2,
          eff_tax2 = 0.2, t_n = 100)

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0,
          C_X = 0, C_R1 = 0, C_R2 =0)

data_store <- data.frame(matrix(ncol = 6, nrow = length(comb_sigma)))

for(i in 1:nrow(comb_sigma)){
  parms1[c("sigma_1", "sigma_2")] <- comb_sigma[i,]
  
  run <- data.frame(ode(y = init, func = amr, times = seq(0, 365), parms = parms1))
  run[9:11] <- run[,3:5]/rowSums(run[,3:5])
  run$avgRes <- rowMeans(run[,10:11])
  colnames(run)[9:12] <- c("WT", "PropRes1", "PropRes2", "AvgPropRes")
  
  Int_cum <- data.frame("time" = run$time,
                        "cum_resprop" = rowMeans(run[,7:8])/2,
                        "cum_sum" = rowSums(run[,6:8]))
  
  print(i)
  
  data_store[i,] <- c(comb_sigma[i,],
                      tail(Int_cum$cum_sum,1),
                      tail(Int_cum$cum_resprop,1),
                      tail(run[,7],1),
                      tail(run[,8],1))
}

colnames(data_store) <- c("sigma1", "sigma2", "wt_res", "avgres", "res1", "res2")

#Plotting Scenarios

plot_list <- list()

for(i in 1:4){
  plot_list[[i]] <- ggplot(data_store, aes_string(x = "sigma1", y = "sigma2", 
                                           fill = c("wt_res", "avgres", "res1", "res2")[i])) + 
    geom_tile() + labs(fill = c("Tot_Inf", "AvgRes", "Res1", "Res 2")[i]) +
    scale_x_discrete(name = bquote("Fraction Pop Treated with Class 1 ("*sigma[1]*")"), expand = c(0, 0)) +  
    scale_y_discrete(name = bquote("Fraction Pop Treated with Class 2 ("*sigma[2]*")"), expand = c(0, 0)) +
    scale_fill_viridis() 
}

ggarrange(plotlist = plot_list)

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

# Parameters and Initial Conditions ---------------------------------------

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0)

parms = c(beta = 0.600094, sigma_1 = 0.35, sigma_2 = 0.35, mu_w = 1.128800, mu_r = 1.247600, 
          mu_t = 1.366400, eta_wr = 0.112880, eta_rw= 0.124760, c1= 0.680032, c2 = 0.560044,
          eff_tax1 = 0.620038, eff_tax2 = 0.680032, t_n = 1000)

# Run the Model -----------------------------------------------------------

testrun <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms))

additions <- length(grep("R",colnames(testrun)))
testrun$X[is.na(testrun$X)] <- 1
testrun[3:(3+additions)][is.na(testrun[3:(3+additions)]) | testrun[3:(3+additions)] < 0 ] <- 0
testrun$AvgRes <- rowMeans(testrun[,4:5])

comb_prop = melt(testrun, id.vars = "time", measure.vars = c("X","W","R1","R2","AvgRes"))

ggplot(data = comb_prop, aes(x = time, y = value, color = variable)) + geom_line() + theme_bw() +
  scale_x_continuous(name = "Time", expand = c(0, 0)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,1)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom")

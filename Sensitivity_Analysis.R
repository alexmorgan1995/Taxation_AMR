library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis")
library("sensitivity")

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

# Functions for the Outcome Measures --------------------------------------

#The first Stage is to create a wrapper function for the eFAST analysis 
#The purpose of this is to have a function which can be used to take in a model input and to output the exact outcome measure you are looking to explore. 

ode_function <- function(x, init, outcome) {
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms = c(beta = x[z,1], sigma_1 = 0.35, sigma_2 = 0.35, mu_x = x[z,2], mu_r1 = x[z,3],
              mu_r2 = x[z,4], mu_t = x[z,5], eta_xr_1 = x[z,6], eta_rx_1 = x[z,7],
              eta_xr_2 = x[z,8], eta_rx_2 = x[z,9], eff_tax1 = x[z,10],
              eff_tax2 = x[z,11], t_n = 500)

    run <- data.frame(ode(y = init, func = amr, times = seq(0, 1000), parms = parms))
    
    print(paste0(round(z/nrow(x), digits  = 3)*100,"%"))
    
    return_vec[z] <- integral(run, parms["t_n"])[outcome]
  }
  return(return_vec)
}

#Test the output of the function
ode_function(data.frame(beta = 0.5, sigma_1 = 0.35, sigma_2 = 0.35, mu_x = 1/15, mu_r1 = 1/10,
               mu_r2 = 1/12, mu_t = 1/4, eta_xr_1 = 0.01, eta_rx_1 = 0.01,
               eta_xr_2 = 0.01, eta_rx_2 = 0.01, eff_tax1 = .5,
               eff_tax2 = .5, t_n = 500), 
             c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0),
              2)

# Run the FAST -----------------------------------------------------------

factors <- c("beta","mu_x", "mu_r1",
             "mu_r2", "mu_t", "eta_xr_1", "eta_rx_1",
             "eta_xr_2", "eta_rx_2", "eff_tax1",
             "eff_tax2")

start_time <- Sys.time()

init <- c(S = 0.99, X = 1-0.99, R1 = 0, R2 =0)

testfbd <- fast99(model = ode_function, factors = factors, n = 2500, 
                  q.arg = list(list(min=0.05, max=5), 
                               list(min=1/100, max=1),
                               list(min=1/100, max = 1),
                               list(min=1/100, max=1),
                               list(min=1/100, max=1),
                               list(min=0.001, max=0.1),
                               list(min=0.001, max=0.1),
                               list(min=0.001, max=0.1),
                               list(min=0.001, max=0.1),
                               list(min=0, max=1),
                               list(min=0, max=1)),
                  init = init,
                  outcome = 2)

end_time <- Sys.time()
end_time - start_time

# Plotting the Output -----------------------------------------------------

par(mfrow = c(2,1))

plot(testfbd)

f99_dataframe_fbd <- data.frame(X=colnames(testfbd$X), testfbd$D1/testfbd$V, 1 - testfbd$Dt/testfbd$V - testfbd$D1/testfbd$V , 1 - testfbd$Dt/testfbd$V)
colnames(f99_dataframe_fbd)[-1] <- c("first.order","higher.order", "total.order")

plotdata_FBD <- melt(f99_dataframe_fbd, id.vars = "X", measure.vars = c("higher.order","first.order"))

plotdata_FBD$X <- factor(plotdata_FBD$X, levels = reorder(unique(plotdata_FBD$X), - f99_dataframe_fbd$total.order))

# Plotting the eFAST ------------------------------------------------------

p_efast_fbd <- ggplot(plotdata_FBD, aes(fill = variable, x =  X, y = value)) + theme_bw() + 
  geom_col(color = "black",position= "stack") + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(legend.position=c(0.25, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5)) + 
  scale_fill_manual(labels = c("Second Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = bquote(bold("eFAST total/first order sensitivity indices")),
       x ="", y = "Sensitivity Index")

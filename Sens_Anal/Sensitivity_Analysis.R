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

# Functions for the Outcome Measures --------------------------------------

#The first Stage is to create a wrapper function for the eFAST analysis 
#The purpose of this is to have a function which can be used to take in a model input and to output the exact outcome measure you are looking to explore. 

ode_function <- function(x, init, outcome) {
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms = c(beta = x[z,1], sigma_1 = 0.35, sigma_2 = 0.35, mu_w = x[z,2],
              mu_r = x[z,3], mu_t = x[z,4], eta_wr = x[z,5], eta_rw = x[z,6],
              c1= x[z,7], c2 = x[z,8], eff_tax1 = x[z,9], eff_tax2 = x[z,10], 
              t_n = 1000)

    run <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms))
    
    #Cleaning the ODE outpout to ensure that there are no NAs or values below 
    
    additions <- length(grep("R",colnames(run)))
    run$X[is.na(run$X)] <- 1
    run[3:(3+additions)][is.na(run[3:(3+additions)]) | run[3:(3+additions)] < 0 ] <- 0
    run$AvgRes <- rowMeans(run[,4:5])
    
    #Progress Bar
    print(paste0(round(z/nrow(x), digits  = 4)*100,"%"))
    
    return_vec[z] <- integral(run, parms["t_n"])[outcome]
    if(is.na(integral(run, parms["t_n"])[outcome])) {
      print(parms)
    }
  }
  return(return_vec)
}

#Test the output of the function
ode_function(data.frame(beta = 2, mu_w = 1/15, mu_r = 1/10, 
                        mu_t = 1/5, eta_wr = 0.02, eta_rw= 0.02, c1= 0.9, c2 = 0.85,
                        eff_tax1 = 1, eff_tax2 = 1), 
             c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0),
              1)

# Run the FAST -----------------------------------------------------------

factors <- c("beta","mu_w", "mu_r", "mu_t", "eta_wr", "eta_rw",
             "c1", "c2", "eff_tax1", "eff_tax2")

start_time <- Sys.time()

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0)

testfbd <- fast99(model = ode_function, factors = factors, n = 200, 
                  q.arg = list(list(min=0.0001, max=10), 
                               list(min=1/50, max=1/0.5),
                               list(min=1/50, max = 1/0.5),
                               list(min=1/50, max=1/0.5),
                               list(min=0.002, max=0.2),
                               list(min=0.002, max=0.2),
                               list(min=0.0001, max=1),
                               list(min=0.0001, max=1),
                               list(min=0.0001, max=1),
                               list(min=0.0001, max=1)),
                  init = init,
                  outcome = 2)


testfbd$y
end_time <- Sys.time()
end_time - start_time

# Plotting the Output -----------------------------------------------------

par(mfrow = c(2,1))
summary(testfbd[1])
plot(testfbd, ylim= c(0, 3))

testfbd$y

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



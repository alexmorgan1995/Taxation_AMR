library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  data_temp <- data_temp[,-14]
  data_temp$Res1 <- rowSums(data_temp[,5:6])
  data_temp$Res2 <- rowSums(data_temp[,7:8])
  data_temp$Res3 <- rowSums(data_temp[,9:10])

  out_vec <- signif(c(sum(data_temp[3:10]),
                      sum(rowMeans(data_temp[14:16]))),
                    5)
  return(out_vec)
}

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    R_base1 <- R1
    R_base2 <- R2
    
    if(t > t_n) {
      R1 <- R_base1*(1-eff_tax1_1) 
      R2 <- R_base2*(1-eff_tax2_1)
    }
    
    if(t > (t_n + time_between)) {
      R1 <- R_base1*(1-eff_tax1_2) 
      R2 <- R_base2*(1-eff_tax2_2)
    }
    
    if(t > (t_n + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_3) 
      R2 <- R_base2*(1-eff_tax2_3)
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_4) 
      R2 <- R_base2*(1-eff_tax2_4)
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_5) 
      R2 <- R_base2*(1-eff_tax2_5)
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_6) 
      R2 <- R_base2*(1-eff_tax2_6)
    }
    
    R1 <- ifelse(R1 > 0, R1, 0)
    R2 <- ifelse(R2 > 0, R2, 0)
    
  dS = lambda - lambda*S + 
    mu_wt*eta*Wc + mu_r*eta*Rc1 + mu_r*eta*Rc2 + 
    P1*theta + P2*theta + 
    Pr*rho*(R1+R2)*Wc + 
    Pr*(R2)*Rc1 + Pr*(R1)*Rc2 - Pr*(R1+R2)*S - 
    beta*S*(Wc + Wi) - beta*S*(Ri1 + Rc1)*fc1 - beta*S*(Ri2 + Rc2)*fc2 +
    tau*(R1+R2)*Wi + tau*(R2)*Ri1 + tau*(R1)*Ri2
    
  dWc = -lambda*Wc - mu_wt*eta*Wc + mu_wt*Wi - Pr*(R1+R2)*zeta*Wc + 
    beta*S*(Wc + Wi)*psi + kappa*(Rc1 + Rc2) 
  
  dWi = -lambda*Wi- mu_wt*Wi - tau*(R1+R2)*Wi + beta*S*(Wc + Wi)*(1-psi)
  
  dRc1 = -lambda*Rc1 - mu_r*eta*Rc1 + mu_r*Ri1 + Pr*R1*(1-zeta)*Wc - Pr*(R2)*zeta*Rc1 + 
    beta*S*(Ri1 + Rc1)*psi*fc1 - kappa*Rc1 + Pr*R1*(1-rho)*S
  
  dRi1 = -lambda*Ri1 - mu_r*Ri1 - tau*(R2)*Ri1 + beta*S*(Ri1 + Rc1)*(1-psi)*fc1
  
  dRc2 = -lambda*Rc2 - mu_r*eta*Rc2 + mu_r*Ri2 + Pr*R2*(1-zeta)*Wc - Pr*(R1)*zeta*Rc2 + 
    beta*S*(Ri2 + Rc2)*psi*fc2 - kappa*Rc2 + Pr*R2*(1-rho)*S
  
  dRi2 = -lambda*Ri2 - mu_r*Ri2 - tau*(R1)*Ri2 + beta*S*(Ri2 + Rc2)*(1-psi)*fc2
  
  dP1 = -lambda*P1 + Pr*R1*rho*S - P1*theta
  
  dP2 = -lambda*P2 + Pr*R2*rho*S - P2*theta
  
    #Calculating the Proportion Integrals
    
    return(list(c(dS, dWc, dWi, dRc1, dRi1, dRc2, dRi2, dP1, dP2)))
  })
}

# Run the ODE -------------------------------------------------------------

#Init Parameters 
init <- c(S = 0.99, Wc = 1-0.99, Wi = 0, 
          Rc1 = 0, Ri1 = 0,
          Rc2 = 0, Ri2 = 0,
          P1 = 0, P2 = 0)

parms = c(beta = 5, R1 = 0.3, R2 = 0.3, 
          lambda = 0.1, Pr = 10, rho = 0.01, theta = 0.01,
          tau = 0.01, psi = 0.95, kappa = 0.1, eta = 0.01, 
          fc1 = 0.9, fc2 = 0.86, 
          mu_wt = 1/12, mu_r = 1/10,  
          zeta = 0.1, 
          
          eff_tax1_1 = 0, eff_tax2_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0,
          
          t_n = 1000, time_between = Inf)


#Flat Tax
parms1 <- parms
parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))

testrun_flat$WT <- rowSums(testrun_flat[,3:4])
testrun_flat$Res1 <- rowSums(testrun_flat[,5:6])
testrun_flat$Res2 <- rowSums(testrun_flat[,7:8])
testrun_flat$proph <- rowSums(testrun_flat[,9:10])

t_melt <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[11:14])
#t_melt <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[-1])

ggplot(t_melt, aes(time, value, color = variable)) + geom_line(size = 1.2) + theme_bw() + 
  scale_x_continuous(name = "Time (days)", expand = c(0, 0)) +  scale_y_continuous(name = "Proportion Infected", expand = c(0, 0), limits = c(0,1)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") 

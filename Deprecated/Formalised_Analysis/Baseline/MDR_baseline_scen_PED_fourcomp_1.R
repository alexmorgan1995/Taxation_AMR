library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
library("rbenchmark")

rm(list=ls())

# ODEs --------------------------------------------------------------------


amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    sigma_base1 <- sigma1; sigma_base2 <- sigma2; sigma_base3 <- sigma3; sigma_base4 <- sigma4
    
    if(t > t_n) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,1]*PED[1,1]) + (eff_tax[2,1]*PED[2,1]) + (eff_tax[3,1]*PED[3,1]) + (eff_tax[4,1]*PED[4,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,1]*PED[1,2]) + (eff_tax[2,1]*PED[2,2]) + (eff_tax[3,1]*PED[3,2]) + (eff_tax[4,1]*PED[4,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,1]*PED[1,3]) + (eff_tax[2,1]*PED[2,3]) + (eff_tax[3,1]*PED[3,3]) + (eff_tax[4,1]*PED[4,3])))
      sigma4 <- sigma_base4*(1 + ((eff_tax[1,1]*PED[1,4]) + (eff_tax[2,1]*PED[2,4]) + (eff_tax[3,1]*PED[3,4]) + (eff_tax[4,1]*PED[4,4])))
    }
    
    if(t > (t_n + time_between)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,2]*PED[1,1]) + (eff_tax[2,2]*PED[2,1]) + (eff_tax[3,2]*PED[3,1]) + (eff_tax[4,2]*PED[4,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,2]*PED[1,2]) + (eff_tax[2,2]*PED[2,2]) + (eff_tax[3,2]*PED[3,2]) + (eff_tax[4,2]*PED[4,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,2]*PED[1,3]) + (eff_tax[2,2]*PED[2,3]) + (eff_tax[3,2]*PED[3,3]) + (eff_tax[4,2]*PED[4,3])))
      sigma4 <- sigma_base4*(1 + ((eff_tax[1,2]*PED[1,4]) + (eff_tax[2,2]*PED[2,4]) + (eff_tax[3,2]*PED[3,4]) + (eff_tax[4,2]*PED[4,4])))
    }
    
    if(t > (t_n + time_between*2)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,3]*PED[1,1]) + (eff_tax[2,3]*PED[2,1]) + (eff_tax[3,3]*PED[3,1]) + (eff_tax[4,3]*PED[4,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,3]*PED[1,2]) + (eff_tax[2,3]*PED[2,2]) + (eff_tax[3,3]*PED[3,2]) + (eff_tax[4,3]*PED[4,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,3]*PED[1,3]) + (eff_tax[2,3]*PED[2,3]) + (eff_tax[3,3]*PED[3,3]) + (eff_tax[4,3]*PED[4,3])))
      sigma4 <- sigma_base4*(1 + ((eff_tax[1,3]*PED[1,4]) + (eff_tax[2,3]*PED[2,4]) + (eff_tax[3,3]*PED[3,4]) + (eff_tax[4,3]*PED[4,4])))
    }
    
    if(t > (t_n + time_between*3)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,4]*PED[1,1]) + (eff_tax[2,4]*PED[2,1]) + (eff_tax[3,4]*PED[3,1]) + (eff_tax[4,4]*PED[4,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,4]*PED[1,2]) + (eff_tax[2,4]*PED[2,2]) + (eff_tax[3,4]*PED[3,2]) + (eff_tax[4,4]*PED[4,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,4]*PED[1,3]) + (eff_tax[2,4]*PED[2,3]) + (eff_tax[3,4]*PED[3,3]) + (eff_tax[4,4]*PED[4,3])))
      sigma4 <- sigma_base4*(1 + ((eff_tax[1,4]*PED[1,4]) + (eff_tax[2,4]*PED[2,4]) + (eff_tax[3,4]*PED[3,4]) + (eff_tax[4,4]*PED[4,4])))
    }
    
    if(t > (t_n + time_between*4)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,5]*PED[1,1]) + (eff_tax[2,5]*PED[2,1]) + (eff_tax[3,5]*PED[3,1]) + (eff_tax[4,5]*PED[4,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,5]*PED[1,2]) + (eff_tax[2,5]*PED[2,2]) + (eff_tax[3,5]*PED[3,2]) + (eff_tax[4,5]*PED[4,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,5]*PED[1,3]) + (eff_tax[2,5]*PED[2,3]) + (eff_tax[3,5]*PED[3,3]) + (eff_tax[4,5]*PED[4,3])))
      sigma4 <- sigma_base4*(1 + ((eff_tax[1,5]*PED[1,4]) + (eff_tax[2,5]*PED[2,4]) + (eff_tax[3,5]*PED[3,4]) + (eff_tax[4,5]*PED[4,4])))
    }
    
    if(t > (t_n + time_between*5)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,6]*PED[1,1]) + (eff_tax[2,6]*PED[2,1]) + (eff_tax[3,6]*PED[3,1]) + (eff_tax[4,6]*PED[4,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,6]*PED[1,2]) + (eff_tax[2,6]*PED[2,2]) + (eff_tax[3,6]*PED[3,2]) + (eff_tax[4,6]*PED[4,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,6]*PED[1,3]) + (eff_tax[2,6]*PED[2,3]) + (eff_tax[3,6]*PED[3,3]) + (eff_tax[4,6]*PED[4,3])))
      sigma4 <- sigma_base4*(1 + ((eff_tax[1,6]*PED[1,4]) + (eff_tax[2,6]*PED[2,4]) + (eff_tax[3,6]*PED[3,4]) + (eff_tax[4,6]*PED[4,4])))
    }
    
    sigma1 <- ifelse(sigma1 > 0, sigma1, 0)
    sigma2 <- ifelse(sigma2 > 0, sigma2, 0)
    sigma3 <- ifelse(sigma3 > 0, sigma3, 0)    
    sigma4 <- ifelse(sigma4 > 0, sigma4, 0)
    
    if((sigma1 + sigma2 + sigma3 + sigma4) > 1) {
      sigma1 <- sigma1/(sigma1 + sigma2 + sigma3 + sigma4)
      sigma2 <- sigma2/(sigma1 + sigma2 + sigma3 + sigma4)
      sigma3 <- sigma3/(sigma1 + sigma2 + sigma3 + sigma4)
      sigma4 <- sigma4/(sigma1 + sigma2 + sigma3 + sigma4)
    }
    
    dX = lambda - lambda*X - 
      beta*X*(Wt + R1*c1 + R2*c2 + R3*c3 + R4*c4 + 
                R12*c12 + R13*c13 + R14*c14 + 
                R23*c23 + R24*c24 + R34*c34 + 
                R123*c123 + R124*c124 + R134*c134 + R234*c234 +
                R1234*c1234) +
      r_wt*Wt*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + 
      r_r*R1*sigma1 + r_r*R2*sigma2 + r_r*R3*sigma3 + r_r*R4*sigma4 +
      r_rr*R12*(sigma1 + sigma2) + r_rr*R13*(sigma1 + sigma3) + r_rr*R23*(sigma2 + sigma3) + r_rr*R14*(sigma1 + sigma4) + 
      r_rr*R24*(sigma2 + sigma4) + r_rr*R34*(sigma3 + sigma4) + 
      r_rrr*R123*(sigma1 + sigma2 + sigma3) + r_rrr*R124*(sigma1 + sigma2 + sigma4) + 
      r_rrr*R234*(sigma2 + sigma3 + sigma4) + r_rrr*R134*(sigma1 + sigma3 + sigma4) + 
      r_rrrr*R1234*(sigma1 + sigma2 + sigma3 + sigma4) + 
      r_t*(1-rho)*(Wt*(sigma1 + sigma2 + sigma3 + sigma4) + 
                     R1*(sigma2 + sigma3 + sigma4) + 
                     R2*(sigma1 + sigma3 + sigma4) + 
                     R3*(sigma1 + sigma2 + sigma4) + 
                     R4*(sigma1 + sigma2 + sigma3) + 
                     R12*(sigma3 + sigma4) + R13*(sigma2 + sigma4) + 
                     R14*(sigma2 + sigma3) + R23*(sigma1 + sigma4) +
                     R24*(sigma1 + sigma3) + R34*(sigma1 + sigma2) +
                     R123*(sigma4) + R124*(sigma3) + 
                     R234*(sigma1) + R134*(sigma2))
    
    dWt = - lambda*Wt + beta*X*Wt - 
      r_wt*Wt*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) - 
      r_t*Wt*(1-rho)*(sigma1 + sigma2 + sigma3 + sigma4) +
      eta_rw*(R1 + R2 + R3 + R4 + 
                R12 + R13 + R14 + R23 + R24 + R34 +
                R123 + R124 + R134 + R234 + 
                R1234)*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) - 
      eta_wr*Wt*rho*(sigma1 + sigma2 + sigma3 + sigma4)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma2 + sigma3 + sigma4)*R1 - r_r*sigma1*R1 - eta_rr*R1*rho*sigma2 -
      eta_rr*R1*rho*sigma3 - eta_rr*R1*rho*sigma4 - eta_rw*R1*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_wr*rho*Wt*sigma1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma1 + sigma3 + sigma4)*R2 - r_r*sigma2*R2 - eta_rr*R2*rho*sigma1 -
      eta_rr*R2*rho*sigma3 - eta_rr*R2*rho*sigma4 - eta_rw*R2*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_wr*rho*Wt*sigma2
    
    dR3 = - lambda*R3 + beta*X*R3*c3 - r_t*(1-rho)*(sigma1 + sigma2 + sigma4)*R3 - r_r*sigma3*R3 - eta_rr*R3*rho*sigma1 -
      eta_rr*R3*rho*sigma2 - eta_rr*R3*rho*sigma4 - eta_rw*R3*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_wr*rho*Wt*sigma3
    
    dR4 = - lambda*R4 + beta*X*R4*c4 - r_t*(1-rho)*(sigma1 + sigma2 + sigma3)*R4 - r_r*sigma4*R4 - eta_rr*R4*rho*(sigma1 + sigma2 + sigma3) - 
      eta_rw*R4*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_wr*rho*Wt*sigma4
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_t*(1-rho)*(sigma3 + sigma4)*R12 - r_rr*(sigma1 + sigma2)*R12 - eta_rrr*R12*rho*(sigma3) - eta_rrr*R12*rho*sigma4 -
      eta_rw*R12*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rr*R1*rho*sigma2 + eta_rr*R2*rho*sigma1 
    
    dR13 = - lambda*R13 + beta*X*R13*c13 - r_t*(1-rho)*(sigma2 + sigma4)*R13 - r_rr*(sigma1 + sigma3)*R13 - eta_rrr*R13*rho*(sigma2) - eta_rrr*R13*rho*sigma4 -
      eta_rw*R13*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rr*R1*rho*sigma3 + eta_rr*R3*rho*sigma1 
    
    dR14 = - lambda*R14 + beta*X*R14*c14 - r_t*(1-rho)*(sigma2 + sigma3)*R14 - r_rr*(sigma1 + sigma4)*R14 - eta_rrr*R14*rho*(sigma2) - eta_rrr*R14*rho*sigma3 -
      eta_rw*R14*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rr*R1*rho*sigma4 + eta_rr*R4*rho*sigma1 
    
    dR23 = - lambda*R23 + beta*X*R23*c23 - r_t*(1-rho)*(sigma1 + sigma4)*R23 - r_rr*(sigma2 + sigma3)*R23 - eta_rrr*R23*rho*(sigma1) - eta_rrr*R23*rho*sigma4 -
      eta_rw*R23*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rr*R2*rho*sigma3 + eta_rr*R3*rho*sigma2
    
    dR24 = - lambda*R24 + beta*X*R24*c24 - r_t*(1-rho)*(sigma1 + sigma3)*R24 - r_rr*(sigma2 + sigma4)*R24 - eta_rrr*R24*rho*(sigma1) - eta_rrr*R24*rho*sigma3 -
      eta_rw*R24*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rr*R2*rho*sigma4 + eta_rr*R4*rho*sigma2 
    
    dR34 = - lambda*R34 + beta*X*R34*c34 - r_t*(1-rho)*(sigma1 + sigma2)*R34 - r_rr*(sigma3 + sigma4)*R34 - eta_rrr*R34*rho*(sigma1) - eta_rrr*R34*rho*sigma2 -
      eta_rw*R34*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rr*R3*rho*sigma4 + eta_rr*R4*rho*sigma3 
    
    
    dR123 = - lambda*R123 + beta*X*R123*c123 - r_t*(1-rho)*R123*(sigma4) - r_rrr*(sigma1 + sigma2 + sigma3)*R123 - 
      eta_rw*R123*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rrr*rho*(sigma3*R12) + eta_rrr*rho*sigma2*R13 + eta_rrr*rho*sigma1*R23 - 
      eta_rrrr*rho*R123*(sigma4)
    
    dR124 = - lambda*R124 + beta*X*R124*c124 - r_t*(1-rho)*R124*(sigma3) - r_rrr*(sigma1 + sigma2 + sigma4)*R124 - 
      eta_rw*R124*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rrr*rho*(sigma4*R12) + eta_rrr*rho*sigma2*R14 + eta_rrr*rho*sigma1*R24 - 
      eta_rrrr*rho*R124*(sigma3)
    
    dR134 = - lambda*R134 + beta*X*R134*c134 - r_t*(1-rho)*R134*(sigma2) - r_rrr*(sigma1 + sigma3 + sigma4)*R134 - 
      eta_rw*R134*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rrr*rho*(sigma4*R13) + eta_rrr*rho*sigma1*R34 + eta_rrr*rho*sigma3*R14 - 
      eta_rrrr*rho*R134*(sigma2)
    
    dR234 = - lambda*R234 + beta*X*R234*c234 - r_t*(1-rho)*R234*(sigma1) - r_rrr*(sigma2 + sigma3 + sigma4)*R234 - 
      eta_rw*R234*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + eta_rrr*rho*(sigma4*R23) + eta_rrr*rho*sigma2*R34 + eta_rrr*rho*sigma3*R24 - 
      eta_rrrr*rho*R234*(sigma1)
    
    dR1234 = - lambda*R1234 + beta*X*R1234*c1234 - r_rrrr*(sigma1 + sigma2 + sigma3 + sigma4)*R1234 - eta_rw*R1234*(1 - (sigma1 + sigma2 + sigma3 + sigma4)) + 
      eta_rrrr*rho*(R123*(sigma4) + R124*(sigma3) + R134*(sigma2) + R234*(sigma1))
    
    return(list(c(dX,dWt,
                  dR1,dR2,dR3,dR4,
                  dR12,dR13,dR14,
                  dR23,dR24,dR34,
                  dR123,dR124,dR134,dR234,
                  dR1234)))
  })
}

# Function to Remove NAs and Round Data -----------------------------------

remNA_func <- function(dataframe){
  n_data <- ncol(dataframe)-1
  
  timing <- t(sapply(1:n_data, function(x)  dataframe[max(which(!is.na(dataframe[,x+1]))),]))
  
  if(timing[1,1] != 10000) {
    for(i in 1:n_data){
      dataframe[seq(timing[[1]]+2,10001),i+1] <- timing[i,i+1]
    }
  }
  dataframe[dataframe < 1e-10] <- 0
  return(dataframe)
}


# Integral Function -------------------------------------------------------

integral <- function(data, t_n, thresh){
  
  #Aggregate the Data into Resistance Classes
  data$aggR1 <- data$R1 + data$R12 + data$R13 + data$R14 + data$R123 + data$R124 + data$R134 + data$R1234
  data$aggR2 <- data$R2 + data$R12 + data$R24 + data$R23 + data$R123 + data$R124 + data$R234 + data$R1234
  data$aggR3 <- data$R3 + data$R13 + data$R23 + data$R34 + data$R123 + data$R134 + data$R234 + data$R1234
  data$aggR4 <- data$R4 + data$R14 + data$R24 + data$R34 + data$R124 + data$R134 + data$R234 + data$R1234
  
  #Subset the Dataframe so that only results after the intervention start
  data_temp <- data[data[,1] > t_n,]
  
  #Determine the X% Thresholds that you want to be under
  thresholds <- unlist(data[t_n-1, 19:22]* thresh)
  under_thresh <- sweep(data[data[,1] > t_n,][,19:22], 2, thresholds)
  
  #Calculate the number of days you are under said threshold
  under_50 <- c(nrow(under_thresh[under_thresh$aggR1 < 0,]), 
                nrow(under_thresh[under_thresh$aggR2 < 0,]), 
                nrow(under_thresh[under_thresh$aggR3 < 0,]),
                nrow(under_thresh[under_thresh$aggR4 < 0,]))
  
  #Find the Sum and make each value proportionate to one another 
  prop_vec <- under_50 / sum(under_50)
  prop_vec <- prop_vec[prop_vec != 0]
  
  prop_vec_avganti <- sum(under_50)/7000
  
  #Output the Optimisation Criteria 
  out_vec <- signif(c(sum(data_temp[3:18]),
                      sum(rowMeans(data_temp[19:22])),
                      prop_vec_avganti,
                      -sum(sapply(1:length(prop_vec), function(x) prop_vec[x]*log(prop_vec[x])))), 5)
  
  return(out_vec)
}


# Function to Aggregate Resistance ----------------------------------------

agg_func <- function(data) {
  agg_data <- data.frame("time" = data$time,
                         "Susc" = data$X,
                         "WT" = data$Wt, 
                         "R1" = data$R1 + data$R12 + data$R13 + data$R14 + data$R123 + data$R124 + data$R134 + data$R1234,
                         "R2" = data$R2 + data$R12 + data$R24 + data$R23 + data$R123 + data$R124 + data$R234 + data$R1234,
                         "R3" = data$R3 + data$R13 + data$R23 + data$R34 + data$R123 + data$R134 + data$R234 + data$R1234,
                         "R4" = data$R4 + data$R14 + data$R24 + data$R34 + data$R124 + data$R134 + data$R234 + data$R1234)
  return(agg_data)
}

# Function to Extract Antibiotic Usage Reductions -------------------------

usage_fun <- function(parms){
  usage = data.frame("time" = seq(0,7000),
                     "PopUsage1" = c(sapply(1:6, function(x) 
                       rep(parms[["sigma1"]] * (1 + ((parms[["eff_tax"]][1,x]*parms[["PED"]][1,1]) + (parms[["eff_tax"]][2,x]*parms[["PED"]][2,1]) + (parms[["eff_tax"]][3,x]*parms[["PED"]][3,1]) + (parms[["eff_tax"]][4,x]*parms[["PED"]][4,1]))), 365*3)),
                       parms[["sigma1"]] * rep((1 + ((parms[["eff_tax"]][1,6]*parms[["PED"]][1,1]) + (parms[["eff_tax"]][2,6]*parms[["PED"]][2,1]) + (parms[["eff_tax"]][3,6]*parms[["PED"]][3,1]) + (parms[["eff_tax"]][4,6]*parms[["PED"]][4,1]))), 
                                               7001-365*3*6)),
                     
                     "PopUsage2" = c(sapply(1:6, function(x) 
                       rep(parms[["sigma2"]] * (1 + ((parms[["eff_tax"]][1,x]*parms[["PED"]][1,2]) + (parms[["eff_tax"]][2,x]*parms[["PED"]][2,2]) + (parms[["eff_tax"]][3,x]*parms[["PED"]][3,2]) + (parms[["eff_tax"]][4,x]*parms[["PED"]][4,2]))), 365*3)),
                       parms[["sigma2"]] * rep((1 + ((parms[["eff_tax"]][1,6]*parms[["PED"]][1,2]) + (parms[["eff_tax"]][2,6]*parms[["PED"]][2,2]) + (parms[["eff_tax"]][3,6]*parms[["PED"]][3,2]) + (parms[["eff_tax"]][4,6]*parms[["PED"]][4,2]))), 
                                               7001-365*3*6)),
                     
                     "PopUsage3" = c(sapply(1:6, function(x) 
                       rep(parms[["sigma3"]] * (1 + ((parms[["eff_tax"]][1,x]*parms[["PED"]][1,3]) + (parms[["eff_tax"]][2,x]*parms[["PED"]][2,3]) + (parms[["eff_tax"]][3,x]*parms[["PED"]][3,3]) + (parms[["eff_tax"]][4,x]*parms[["PED"]][4,3]))), 365*3)),
                       parms[["sigma3"]] * rep((1 + ((parms[["eff_tax"]][1,6]*parms[["PED"]][1,3]) + (parms[["eff_tax"]][2,6]*parms[["PED"]][2,3]) + (parms[["eff_tax"]][3,6]*parms[["PED"]][3,3]) + (parms[["eff_tax"]][4,6]*parms[["PED"]][4,3]))), 
                                               7001-365*3*6)),
                     
                     "PopUsage4" = c(sapply(1:6, function(x) 
                       rep(parms[["sigma4"]] * (1 + ((parms[["eff_tax"]][1,x]*parms[["PED"]][1,4]) + (parms[["eff_tax"]][2,x]*parms[["PED"]][2,4]) + (parms[["eff_tax"]][3,x]*parms[["PED"]][3,4]) + (parms[["eff_tax"]][4,x]*parms[["PED"]][4,4]))), 365*3)),
                       parms[["sigma4"]] * rep((1 + ((parms[["eff_tax"]][1,6]*parms[["PED"]][1,4]) + (parms[["eff_tax"]][2,6]*parms[["PED"]][2,4]) + (parms[["eff_tax"]][3,6]*parms[["PED"]][3,4]) + (parms[["eff_tax"]][4,6]*parms[["PED"]][4,4]))), 
                                               7001-365*3*6)))
  
  usage[usage < 0] <- 0
  usage$totusage = rowSums(usage[2:5])
  usage$reduc_use = parms[["sigma1"]] + parms[["sigma2"]] + parms[["sigma3"]] + parms[["sigma4"]] - usage$totusage
  return(usage)
}

# Dual Model --------------------------------------------------------------


multi_int_fun <- function(int_gen, time_between, parms, init, func, agg_func){
  
  parms["time_between"] <- time_between
  
  #First Run
  run_1rd <- agg_func(remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))))
  
  values_1rd <<- tail(run_1rd, 1)[4:7]

  low_char_1rd <- names(values_1rd)[which.min(values_1rd)]
  high_char_1rd <- names(values_1rd)[which.max(values_1rd)]
  med_char_1rd <- names(values_1rd)[setdiff(1:4, c(which.min(values_1rd), which.max(values_1rd)))]
  med_char_val <- mean(as.numeric(values_1rd[med_char_1rd]))

  parms[["eff_tax"]][as.numeric(substr(low_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[low_char_1rd]/med_char_val)))
  parms[["eff_tax"]][as.numeric(substr(high_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[high_char_1rd]/med_char_val)))
  parms[["eff_tax"]][as.numeric(substr(med_char_1rd[1], 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[med_char_1rd][1]/med_char_val)))
  parms[["eff_tax"]][as.numeric(substr(med_char_1rd[2], 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[med_char_1rd][2]/med_char_val)))
  parms[["eff_tax"]][parms[["eff_tax"]] < 0] <- 0
  
  
  if(int_gen > 1) {
    
    for(i in 1:(int_gen-1)) {
      run <- agg_func(remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*i)), parms = parms))))
      values <- tail(run, 1)[4:7]
      
      if(values[1] == 0 & values[2] == 0 & values[3] == 0) {
        parms[["eff_tax"]][,i] <- 0
      } else {
        
        low_char <- names(values)[which.min(values)]
        high_char <- names(values)[which.max(values)]
        med_char <- names(values)[setdiff(1:4, c(which.min(values), which.max(values)))]
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i+1):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/med_char_val)))
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i+1):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/med_char_val)))
        parms[["eff_tax"]][as.numeric(substr(med_char[1], 2, 2)), c((i+1):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[med_char][1],1)/med_char_val)))
        parms[["eff_tax"]][as.numeric(substr(med_char[2], 2, 2)), c((i+1):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[med_char][2],1)/med_char_val)))
        parms[["eff_tax"]][parms[["eff_tax"]] < 0] <- 0
      }
    }
    
  }
  
  out_run <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
  
  return(list(out_run, parms))
}



# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, 
          R1 = 0, R2 = 0, R3 = 0, R4 = 0,
          R12 = 0, R13 = 0, R14 = 0, R23 = 0, R24 = 0, R34 = 0,
          R123 = 0, R124 = 0, R134 = 0, R234 = 0, 
          R1234 = 0)

parms = list(lambda = 1/365*(2), 
             beta = 5, sigma1 = 0.2, sigma2 = 0.2, sigma3 = 0.2, sigma4 = 0.2,
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, r_rrrr = 1/7, 
             r_t = 1/6, eta_wr = 0.3, eta_rw = 0.04, 
             eta_rr = 0.01, eta_rrr = 0.01, eta_rrrr = 0.01,  
             c1 = 0.945, c2 = 0.935, c3 = 0.91, c4 = 0.8,
             c12 = 0.8, c13 = 0.775,  c14 = 0.75, 
             c23 = 0.75, c24 = 0.725, c34 = 0.7,
             c123 = 0.7, c124 = 0.65, c134 = 0.625, c234 = 0.6,
             c1234 = 0.6,
             PED = matrix(c(-1, 0, 0, 0, 
                            0, -1, 0, 0,
                            0, 0, -1, 0,
                            0, 0, 0, -1), #Be aware of this matrix
                          nrow = 4, ncol = 4, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0), 
                              nrow = 4, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

diff_test_run <- agg_func(multi_int_fun(6, 365*3, parms, init, amr, agg_func)[[1]])
diff_test_run <- melt(diff_test_run, id.vars = "time", measure.vars = colnames(diff_test_run)[4:7])
ggplot(diff_test_run, aes(time, value, color = variable)) + geom_line() + theme_bw()


parms1 <- parms
parms1[["eff_tax"]][4,] <- 0.5
test_run <- agg_func(data.frame(remNA_func(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))))
test_run <- melt(test_run, id.vars = "time", measure.vars = colnames(test_run)[4:7])
ggplot(test_run, aes(time, value, color = variable)) + geom_line() + theme_bw()

# Testing the Outcome Measure ---------------------------------------------

parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5

#Baseline Test Run
testrun_flat <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
integral(testrun_flat, 3000, 0.75) 

#Baseline Differential Taxation Run 
diff_test_run <- multi_int_fun(6, 365*3, parms, init, amr, agg_func)[[1]]
integral(diff_test_run, 3000, 0.75) 

# Baseline Model ----------------------------------------------------------

parms1 <- parms
testrun_flat <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
test_run_agg <- agg_func(testrun_flat) 

test_run_agg$AverageRes <- rowMeans(testrun_flat[,4:7])
test_run_agg$TotInf <- rowSums(testrun_flat[,3:18])

test_plot_flat <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[4:7])
ggplot(test_plot_flat, aes(time, value, color = variable)) + geom_line() + theme_bw()

# Intervention Scenarios --------------------------------------------------

#Flat Tax

parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5
testrun_flat <- list(remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))))

#Single Tax 
single_list <- list()
for(i in 1:4) {
  parms1 <- parms
  parms1[["eff_tax"]][i,] <- 0.5
  single_list[[i]] <- data.frame(remNA_func(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  dat <- multi_int_fun(i, 365*3, parms, init, amr, agg_func)[[1]]
  diff_tax_list[[i]] <- data.frame(dat)
                                  
}

# Obtain the Integrals ----------------------------------------------------

list_scen <- unlist(list(testrun_flat, single_list, diff_tax_list), recursive = FALSE)

opt_crit = data.frame("Flat_All" = integral(list_scen[[1]], 3000, 0.5), 
                      "Single_1" = integral(list_scen[[2]], 3000, 0.5),
                      "Single_2" = integral(list_scen[[3]], 3000, 0.5),
                      "Single_3" = integral(list_scen[[4]], 3000, 0.5),
                      "Single_4" = integral(list_scen[[5]], 3000, 0.5),
                      "Diff_1" = integral(list_scen[[6]], 3000, 0.5),
                      "Diff_2" = integral(list_scen[[7]], 3000, 0.5),
                      "Diff_3" = integral(list_scen[[8]], 3000, 0.5),
                      "Diff_4" = integral(list_scen[[9]], 3000, 0.5),
                      "Diff_5" = integral(list_scen[[10]], 3000, 0.5),
                      "Diff_6" = integral(list_scen[[11]], 3000, 0.5))

rownames(opt_crit) <- c("Tot_Inf", "AvgRes","Avg_Anti","Shannons")

rescale_data <- data.frame(t(apply(opt_crit[c(1:3),], MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
rescale_data["Shannon",] <- unlist(c("Flat_All" = NA, (opt_crit[4,2:11] - min(opt_crit[4,2:11]))/diff(range(opt_crit[4,2:11]))))

# Absolute Plotting Matrix ------------------------------------------------

#All of the plots

OG_melt <- melt(as.matrix(opt_crit), measure.vars = colnames(opt_crit))
rescale_melt <- melt(as.matrix(rescale_data), measure.vars = colnames(rescale_data))

OG_melt$rescale <- rescale_melt[,3]

ggplot(OG_melt, aes(Var2, Var1)) + theme_bw() +
  geom_tile(aes(fill = rescale)) + 
  facet_grid(Var1 ~ ., scales = "free_y") +
  geom_text(aes(label=value), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "Intervention", expand = c(0, 0))  +   
  scale_y_discrete(name = "Outcome Measure", expand = c(0, 0)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="none")

# Creating Relative Integrals ---------------------------------------------

diff_tax_list <- list()

for(i in 1:6) {
  diff_tax_list[[i]] <- multi_int_fun(i, 365*3, parms, init, amr, agg_func)[[2]]
}

#Find The Resistance Lost 
#Baseline

base <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms))
base_totinf <- integral(base, 3000, 0.5)[1]
base_avgres_int <- integral(base, 3000, 0.5)[2]

#Lost Resistance Function 
#Supply a dataframe with time and the effective tax over time for class 1, 2 and 3.
#I need to extract the 1-6 data with the 

#Change in Usage From Baseline - Only for differential taxation models 

#Now I need to do this for all 6 scenarios. 
parms_flat <- parms; parms_flat[["eff_tax"]][,]  <- 0.5

#Create Parameter sets for the different single intervention scenarios

end_vals <- tail(agg_func(base)[4:7], 1)

res_order_vec <- c(names(end_vals)[which.max(end_vals)],
                   names(end_vals)[setdiff(1:4, c(which.min(end_vals), which.max(end_vals)))][1],
                   names(end_vals)[setdiff(1:4, c(which.min(end_vals), which.max(end_vals)))][2],
                   names(end_vals)[which.min(end_vals)])

parms_singleHR <- parms; parms_singleHR[["eff_tax"]][as.numeric(substr(res_order_vec[1], 2, 2)), c(1:6)] <- 0.5
parms_singleMR1 <- parms; parms_singleMR1[["eff_tax"]][as.numeric(substr(res_order_vec[2], 2, 2)), c(1:6)] <- 0.5
parms_singleMR2 <- parms; parms_singleMR2[["eff_tax"]][as.numeric(substr(res_order_vec[3], 2, 2)), c(1:6)] <- 0.5
parms_singleLR <- parms; parms_singleLR[["eff_tax"]][as.numeric(substr(res_order_vec[4], 2, 2)), c(1:6)] <- 0.5

parm_list <- list(parms_flat, parms_singleHR, parms_singleMR1, parms_singleMR2, parms_singleLR)

over_parms <- append(parm_list, diff_tax_list)

list_usage <- list()

for(i in 1:11) {
  list_usage[[i]] <- usage_fun(over_parms[[i]])
}

# Relative Reduction Integral Comparison --------------------------------------------

reduc_usage_vec <- sapply(1:length(list_usage), function(x) sum(list_usage[[x]][,6]))
rel_AvgRes <- base_avgres_int - opt_crit[2,]
rel_totinf <- opt_crit[1,] - base_totinf

comp_totinf <- rel_totinf/reduc_usage_vec
comp_res <- rel_AvgRes/reduc_usage_vec

data.frame.rel <- rbind(comp_totinf, comp_res); rownames(data.frame.rel) <- c("Total Infections", "Average Resistance")
data.frame.rel["Available Antibiotics",] <- opt_crit["Avg_Anti",]
data.frame.rel["Shannon",] <- opt_crit["Shannons",]

rescale_data_rel <- data.frame(t(apply(data.frame.rel, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
rescale_data_rel["Available Antibiotics",] <- rescale_data["Avg_Anti",]
rescale_data_rel["Shannon",] <- rescale_data["Shannon",]

OG_melt_rel <- melt(as.matrix(data.frame.rel), measure.vars = colnames(data.frame.rel))
OG_melt_rel$Var2 <- factor(rep(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR1)", "Single Tax (MR2)", "Single Tax (LR)",
                                 "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)", "Diff Tax (4 Rd)",
                                 "Diff Tax (5 Rd)","Diff Tax (6 Rd)"), each = 4))

OG_melt_rel$Var2  <- factor(OG_melt_rel$Var2 , levels = unique(OG_melt_rel$Var2))

rescale_melt_rel <- melt(as.matrix(rescale_data_rel), measure.vars = colnames(rescale_data_rel))

rescale_melt_rel$Var2 <- factor(rep(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR1)", "Single Tax (MR2)", "Single Tax (LR)",
                                      "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)", "Diff Tax (4 Rd)",
                                      "Diff Tax (5 Rd)","Diff Tax (6 Rd)"), each = 4))

rescale_melt_rel$Var2  <- factor(rescale_melt_rel$Var2 , levels = unique(rescale_melt_rel$Var2))

OG_melt_rel$rescale <- rescale_melt_rel[,3]

OG_melt_rel$value <- round(OG_melt_rel$value, digits = 3)

ggplot(OG_melt_rel, aes(Var2, Var1)) + theme_bw() +
  geom_tile(aes(fill = rescale)) + 
  facet_grid(Var1 ~ ., scales = "free_y") +
  geom_text(aes(label=value), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "Intervention", expand = c(0, 0))  +   
  scale_y_discrete(name = "Outcome Measure", expand = c(0, 0)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="none",
        axis.text.x = element_text(angle = 45, hjust=1))

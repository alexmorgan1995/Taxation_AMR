library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:6]),
                      sum(rowMeans(data_temp[4:6]))),5)
  return(out_vec)
}


# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    delta_base1 <- delta1
    delta_base2 <- delta2
    delta_base3 <- delta3
    
    if(t > t_n) {
      delta1 <- delta_base1*eff_tax1_1 
      delta2 <- delta_base2*eff_tax2_1
      delta3 <- delta_base3*eff_tax3_1
    }
    
    if(t > (t_n + time_between)) {
      delta1 <- delta_base1*eff_tax1_2 
      delta2 <- delta_base2*eff_tax2_2
      delta3 <- delta_base3*eff_tax3_2
    }
    
    if(t > (t_n + time_between + time_between)) {
      delta1 <- delta_base1*eff_tax1_3 
      delta2 <- delta_base2*eff_tax2_3
      delta3 <- delta_base3*eff_tax3_3
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      delta1 <- delta_base1*eff_tax1_4 
      delta2 <- delta_base2*eff_tax2_4
      delta3 <- delta_base3*eff_tax3_4
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      delta1 <- delta_base1*eff_tax1_5 
      delta2 <- delta_base2*eff_tax2_5
      delta3 <- delta_base3*eff_tax3_5
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      delta1 <- delta_base1*eff_tax1_6 
      delta2 <- delta_base2*eff_tax2_6
      delta3 <- delta_base3*eff_tax3_6
    }
    
    dX = (1-m-m1-m2-m3-X)*mu + (delta1*tau + delta2*tau + delta3*tau + gamma)*S + 
      (delta2*tau + delta3*tau + gamma)*R1 + (delta1*tau + delta3*tau + gamma)*R2 +
      (delta1*tau + delta2*tau + gamma)*R3 - 
      beta*X*(S + (1-c1)*R1 + (1-c2)*R2 + (1-c3)*R3)
    
    dS = (m-S)*mu - (delta1*tau + delta2*tau + delta3*tau + gamma)*S + beta*S*X + 
      sigma*beta*(c1*R1 + c2*R2 + c3*R3)*S
    
    
    dR1 = (m1-R1)*mu - (delta2*tau + delta3*tau + gamma)*R1 + beta*(1-c1)*R1*X + 
      sigma*beta*(c1*R1 + (c1-c2)*R2 + (c1-c3)*R3)*R1
    
    dR2 = (m2-R2)*mu - (delta1*tau + delta3*tau + gamma)*R2 + beta*(1-c2)*R2*X + 
      sigma*beta*(c2*R2 + (c2-c1)*R1 + (c2-c3)*R3)*R2
    
    dR3 = (m3-R3)*mu - (delta1*tau + delta2*tau + gamma)*R3 + beta*(1-c3)*R3*X + 
      sigma*beta*(c3*R3 + (c3-c1)*R1 + (c3-c2)*R2)*R3
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX, dS, dR1, dR2, dR3)))
  })
}

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, basetax, time_between, parms, init, func){
  
  parms["time_between"] <- time_between
  
  #First Run
  testrun <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))
  values <- tail(testrun, 1)[4:6]
  
  parms[grep("eff_tax1", names(parms), value = T)] <- 1-(basetax*(values[1]/values[2]))
  parms[grep("eff_tax2", names(parms), value = T)] <- 1-(basetax*(values[2]/values[2]))
  parms[grep("eff_tax3", names(parms), value = T)] <- 1-(basetax*(values[3]/values[2]))
  parms[parms < 0] <- 0
  
  testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  
  #Second Run
  if(int_gen >= 2){
    testrun1 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms))
    values1 <- tail(testrun1, 1)[4:6]
    
    low_char <- names(values1)[which.min(values1)]
    high_char <- names(values1)[which.max(values1)]
    med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]
    
    parms[grep(paste0("eff_tax",substr(low_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[low_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[med_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[high_char],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 3){
    testrun2 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms))
    values2 <- tail(testrun2, 1)[4:6]
    
    low_char1 <- names(values2)[which.min(values2)]
    high_char1 <- names(values2)[which.max(values2)]
    med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
    
    parms[grep(paste0("eff_tax",substr(low_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[low_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[med_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[high_char1],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 4){
    testrun3 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms))
    values3 <- tail(testrun3, 1)[4:6]
    
    low_char2 <- names(values3)[which.min(values3)]
    high_char2 <- names(values3)[which.max(values3)]
    med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]
    
    parms[grep(paste0("eff_tax",substr(low_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[low_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[med_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[high_char2],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 5){
    testrun4 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms))
    values4 <- tail(testrun4, 1)[4:6]
    
    low_char3 <- names(values4)[which.min(values4)]
    high_char3 <- names(values4)[which.max(values4)]
    med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
    
    parms[grep(paste0("eff_tax",substr(low_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[low_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[med_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[high_char3],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 6){
    testrun5 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms))
    values5 <- tail(testrun5, 1)[4:6]
    
    low_char4 <- names(values5)[which.min(values5)]
    high_char4 <- names(values5)[which.max(values5)]
    med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
    
    parms[grep(paste0("eff_tax",substr(low_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[low_char4],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[med_char4],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[high_char4],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  
  return(list(testrun, parms))
}


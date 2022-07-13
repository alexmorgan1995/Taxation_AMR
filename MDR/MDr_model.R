library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  data_temp$aggR1 <- data_temp$R1 + data_temp$R12 + data_temp$R13 + data_temp$R123
  data_temp$aggR2 <- data_temp$R2 + data_temp$R12 + data_temp$R23 + data_temp$R123
  data_temp$aggR3 <- data_temp$R3 + data_temp$R13 + data_temp$R23 + data_temp$R123
  
  out_vec <- signif(c(sum(data_temp[3:10]),
                      sum(rowMeans(data_temp[11:13]))), 5)
  
  return(out_vec)
}

# ODEs --------------------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    sigma_base1 <- sigma1
    sigma_base2 <- sigma2
    sigma_base3 <- sigma3
    
    if(t > t_n) {
      sigma1 <- sigma_base1*(1-(eff_tax1_1*PED1)) 
      sigma2 <- sigma_base2*(1-(eff_tax2_1*PED2))
      sigma3 <- sigma_base3*(1-(eff_tax3_1*PED3))
    }
    
    if(t > (t_n + time_between)) {
      sigma1 <- sigma_base1*(1-(eff_tax1_2*PED1)) 
      sigma2 <- sigma_base2*(1-(eff_tax2_2*PED2)) 
      sigma3 <- sigma_base3*(1-(eff_tax3_2*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between)) {
      sigma1 <- sigma_base1*(1-(eff_tax1_3*PED1)) 
      sigma2 <- sigma_base2*(1-(eff_tax2_3*PED2)) 
      sigma3 <- sigma_base3*(1-(eff_tax3_3*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      sigma1 <- sigma_base1*(1-(eff_tax1_4*PED1)) 
      sigma2 <- sigma_base2*(1-(eff_tax2_4*PED2)) 
      sigma3 <- sigma_base3*(1-(eff_tax3_4*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      sigma1 <- sigma_base1*(1-(eff_tax1_5*PED1)) 
      sigma2 <- sigma_base2*(1-(eff_tax2_5*PED2)) 
      sigma3 <- sigma_base3*(1-(eff_tax3_5*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      sigma1 <- sigma_base1*(1-(eff_tax1_6*PED1)) 
      sigma2 <- sigma_base2*(1-(eff_tax2_6*PED2)) 
      sigma3 <- sigma_base3*(1-(eff_tax3_6*PED3)) 
    }
    
    sigma1 <- ifelse(sigma1 > 0, sigma1, 0)
    sigma2 <- ifelse(sigma2 > 0, sigma2, 0)
    sigma3 <- ifelse(sigma3 > 0, sigma3, 0)
    
    
    
    
    dX = lambda - lambda*X - beta*X*(Wt + R1*c1 + R2*c2 + R3*c3 + R12*c12 + R13*c13 + R23*c23 + R123*c123) +
      r_wt*Wt*(1 - sigma1 + sigma2 + sigma3) + r_r*R1*sigma1 + r_r*R2*sigma2 + r_r*R3*sigma3 +
      r_rr*R12*(sigma1 + sigma2) + r_rr*R13*(sigma1 + sigma3) + r_rr*R23*(sigma2 + sigma3) + 
      r_rrr*R123*(sigma1 + sigma2 + sigma3) + r_t*(1-rho)*(Wt*(sigma1 + sigma2 + sigma3) + R1*(sigma2 + sigma3) + 
                                                             R2*(sigma1 + sigma3) + R3*(sigma1 + sigma2) + 
                                                             R12*sigma3 + R13*sigma2 + R23*sigma1)
      
    dWt = - lambda*Wt + beta*X*Wt - r_wt*Wt*(1 - sigma1 + sigma2 + sigma3) - r_t*Wt*(1-rho)*(sigma1 + sigma2 + sigma3) +
      eta_rw*(R1 + R2 + R3 + R12 + R13 + R23 + R123)*(1 - sigma1 + sigma2 + sigma3) - 
      eta_wr*Wt*rho*(sigma1 + sigma2 + sigma3)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma2 + sigma3)*R1 - r_r*sigma1*R1 - eta_rr*R1*rho*sigma2 -
      eta_rr*R1*rho*sigma3 - eta_rw*R1*(1 - sigma1 + sigma2 + sigma3) + eta_wr*rho*Wt*sigma1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma1 + sigma3)*R2 - r_r*sigma2*R2 - eta_rr*R2*rho*sigma1 -
      eta_rr*R2*rho*sigma3 - eta_rw*R2*(1 - sigma1 + sigma2 + sigma3) + eta_wr*rho*Wt*sigma2
    
    dR3 = - lambda*R3 + beta*X*R3*c3 - r_t*(1-rho)*(sigma1 + sigma2)*R3 - r_r*sigma3*R3 - eta_rr*R3*rho*sigma1 -
      eta_rr*R3*rho*sigma2 - eta_rw*R3*(1 - sigma1 + sigma2 + sigma3) + eta_wr*rho*Wt*sigma3
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_t*(1-rho)*sigma3*R12 - r_rr*(sigma1 + sigma2)*R12 - eta_rrr*R12*rho*sigma3 -
      eta_rw*R12*(1 - sigma1 + sigma2 + sigma3) + eta_rr*R1*rho*sigma2 + eta_rr*R2*rho*sigma1 
    
    dR13 = - lambda*R13 + beta*X*R13*c13 - r_t*(1-rho)*sigma2*R13 - r_rr*(sigma1 + sigma3)*R13 - eta_rrr*R13*rho*sigma2 -
      eta_rw*R13*(1 - sigma1 + sigma2 + sigma3) + eta_rr*R1*rho*sigma3 + eta_rr*R3*rho*sigma1 
    
    dR23 = - lambda*R23 + beta*X*R23*c23 - r_t*(1-rho)*sigma1*R23 - r_rr*(sigma2 + sigma3)*R23 - eta_rrr*R23*rho*sigma1 -
      eta_rw*R23*(1 - sigma1 + sigma2 + sigma3) + eta_rr*R2*rho*sigma3 + eta_rr*R3*rho*sigma2 
    
    dR123 = - lambda*R123 + beta*X*R123*c123 - r_rrr*(sigma1 + sigma2 + sigma3)*R123 - eta_rw*R123*(1 - sigma1 + sigma2 + sigma3) + 
      eta_rrr*rho*(sigma3*R12 + sigma2*R13 + sigma1*R23)
    
    return(list(c(dX,dWt,
                  dR1,dR2,dR3,
                  dR12,dR13,dR23,
                  dR123)))
  })
}

# Clean Model Run --------------------------------------------------------------

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

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, time_between, parms, init, func, agg_func){
  
  parms["time_between"] <- time_between
  
  #First Run
  testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms, hmax = 1)))
  testrun <- agg_func(testrun)
  values <- tail(testrun, 1)[4:6]

  parms[grep("eff_tax1", names(parms), value = T)] <- (parms[["base_tax"]]*(values[1]/values[2]))
  parms[grep("eff_tax2", names(parms), value = T)] <- (parms[["base_tax"]]*(values[2]/values[2]))
  parms[grep("eff_tax3", names(parms), value = T)] <- (parms[["base_tax"]]*(values[3]/values[2]))
  parms[parms < 0] <- 0
  
  if(int_gen <= 2 | int_gen >= 6) {
    testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
  }
  
  #Second Run
  if(int_gen >= 2){
    
    testrun1 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms, hmax = 1)))
    testrun1 <- agg_func(testrun1)
    values1 <- tail(testrun1, 1)[4:6]

    if(values1[1] == 0 & values1[2] == 0 & values1[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,3)]] <- 0
    } else {
      low_char <- names(values1)[which.min(values1)]
      high_char <- names(values1)[which.max(values1)]
      med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]
      parms[grep(paste0("eff_tax", substr(low_char, 2, 2)), names(parms), value = T)[-c(1)]] <- (parms[["base_tax"]]*(tail(testrun1[low_char],1)/values[2]))
      parms[grep(paste0("eff_tax", substr(med_char, 2, 2)), names(parms), value = T)[-c(1)]] <- (parms[["base_tax"]]*(tail(testrun1[med_char],1)/values[2]))
      parms[grep(paste0("eff_tax", substr(high_char, 2, 2)), names(parms), value = T)[-c(1)]] <- (parms[["base_tax"]]*(tail(testrun1[high_char],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 2) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 3){
    testrun2 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms, hmax = 1)))
    testrun2 <- agg_func(testrun2)
    values2 <- tail(testrun2, 1)[4:6]
    
    if(values2[1] == 0 & values2[2] == 0 & values2[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,6)]] <- 0
    } else {
      low_char1 <- names(values2)[which.min(values2)]
      high_char1 <- names(values2)[which.max(values2)]
      med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
      parms[grep(paste0("eff_tax",substr(low_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (parms[["base_tax"]]*(tail(testrun2[low_char1],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (parms[["base_tax"]]*(tail(testrun2[med_char1],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (parms[["base_tax"]]*(tail(testrun2[high_char1],1)/values[2]))
      
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 3) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 4){
    testrun3 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms, hmax = 1)))
    testrun3 <- agg_func(testrun3)
    values3 <- tail(testrun3, 1)[4:6]
    
    if(values3[1] == 0 & values3[2] == 0 & values3[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,9)]] <- 0
    } else {
      low_char2 <- names(values3)[which.min(values3)]
      high_char2 <- names(values3)[which.max(values3)]
      med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]
      parms[grep(paste0("eff_tax",substr(low_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (parms[["base_tax"]]*(tail(testrun3[low_char2],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (parms[["base_tax"]]*(tail(testrun3[med_char2],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (parms[["base_tax"]]*(tail(testrun3[high_char2],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 4) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 5){
    testrun4 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms, hmax = 1)))
    testrun4 <- agg_func(testrun4)
    values4 <- tail(testrun4, 1)[4:6]
    
    if(values4[1] == 0 & values4[2] == 0 & values4[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,12)]] <- 0
    } else {
      low_char3 <- names(values4)[which.min(values4)]
      high_char3 <- names(values4)[which.max(values4)]
      med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
      parms[grep(paste0("eff_tax",substr(low_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (parms[["base_tax"]]*(tail(testrun4[low_char3],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (parms[["base_tax"]]*(tail(testrun4[med_char3],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (parms[["base_tax"]]*(tail(testrun4[high_char3],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 5) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 6){
    testrun5 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms, hmax = 1)))
    testrun5 <- agg_func(testrun5)
     values5 <- tail(testrun5, 1)[4:6]
    
    if(values5[1] == 0 & values5[2] == 0 & values5[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,15)]] <- 0
    } else {
      low_char4 <- names(values5)[which.min(values5)]
      high_char4 <- names(values5)[which.max(values5)]
      med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
      parms[grep(paste0("eff_tax",substr(low_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (parms[["base_tax"]]*(tail(testrun5[low_char4],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (parms[["base_tax"]]*(tail(testrun5[med_char4],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (parms[["base_tax"]]*(tail(testrun5[high_char4],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 6) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  return(list(testrun, parms))
}

# Aggregated Function -----------------------------------------------------

agg_func <- function(data) {
  agg_data <- data.frame("time" = data$time,
                         "Susc" = data$X,
                         "WT" = data$Wt, 
                         "R1" = data$R1 + data$R12 + data$R13 + data$R123,
                         "R2" = data$R2 + data$R12 + data$R23 + data$R123,
                         "R3" = data$R3 + data$R13 + data$R23 + data$R123)
  return(agg_data)
}

# Extract Usage -----------------------------------------------------------

usage_fun <- function(parms){
  usage = data.frame("time" = seq(0,7000),
                     "PopUsage1" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax1_", x)]])*parms[["sigma1"]], 365*3)),
                                     rep((1-parms[["eff_tax1_6"]])*parms[["sigma1"]], 7001-365*3*6)),
                     "PopUsage2" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax2_", x)]])*parms[["sigma2"]], 365*3)),
                                     rep((1-parms[["eff_tax2_6"]])*parms[["sigma2"]], 7001-365*3*6)),
                     "PopUsage3" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax3_", x)]])*parms[["sigma3"]], 365*3)),
                                     rep((1-parms[["eff_tax3_6"]])*parms[["sigma3"]], 7001-365*3*6)))
  usage[usage < 0] <- 0
  usage$totusage = rowSums(usage[2:4])
  usage$reduc_use = parms[["sigma1"]] + parms[["sigma2"]] + parms[["sigma3"]] - usage$totusage
  return(usage)
}


# Parameters --------------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = c(lambda = 1/365*(2), 
          beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
          r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
          r_t = 1/7, eta_wr = 0.2, eta_rw = 0.04, 
          eta_rr = 0.1, eta_rrr = 0.1,  
          c1 = 0.95, c2 = 0.92, c3 = 0.85,
          c12 = 0.85, c13 = 0.82, c23 = 0.75,
          c123 = 0.7,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 3000, time_between = Inf, rho = 0.1, base_tax = 0.5)

# Run Model ---------------------------------------------------------------



parms = c(lambda = 1/365*(2), 
          beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
          r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
          r_t = 1/7, eta_wr = 0.5, eta_rw = 0.03, 
          eta_rr = 0.1, eta_rrr = 0.1,  
          c1 = 0.95, c2 = 0.92, c3 = 0.85,
          c12 = 0.85, c13 = 0.82, c23 = 0.75,
          c123 = 0.7,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 3000, time_between = Inf, rho = 0.1, base_tax = 0.5)


parms[grep("eff_tax3", names(parms),value = T)] <- 0.5

testrun <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms)))

int_test <- integral(testrun, 3000)

# Test Plotting -----------------------------------------------------------

test_plot <- melt(testrun, id.vars = "time", measure.vars = colnames(testrun)[-1])

# Aggregated Test Plotting ------------------------------------------------

agg_data <- agg_func(testrun)
agg_data$avgres <- rowMeans(agg_data[,4:6])

test_plot_agg <- melt(agg_data, id.vars = "time", measure.vars = colnames(agg_data)[-1])

ggplot(test_plot_agg, aes(time, value, color = variable)) + geom_line(size = 1.1) + theme_bw() + 
  scale_y_continuous(name = "Average Resistance", expand = c(0, 0)) + scale_x_continuous(name = "Time (days)", expand = c(0, 0)) +
  labs(color="Strains")

# Intervention Scenarios --------------------------------------------------

parms = c(lambda = 1/365*(2), 
          beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
          r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
          r_t = 1/7, eta_wr = 0.5, eta_rw = 0.03, 
          eta_rr = 0.1, eta_rrr = 0.1,  
          c1 = 0.95, c2 = 0.92, c3 = 0.85,
          c12 = 0.85, c13 = 0.82, c23 = 0.75,
          c123 = 0.7,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 3000, time_between = Inf, rho = 0.1, base_tax = 0.5)

#Flat Tax
parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))

#Single Tax 
single_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", i), names(parms1), value = TRUE)] <- 0.5
  single_list[[i]] <- data.frame(remNA_func(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  parms1 <- parms
  dat <- multi_int_fun(i, 365*3, parms1, init, amr, agg_func)[[1]]
  diff_tax_list[[i]] <- dat
}

# Obtain the Integrals ----------------------------------------------------

parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5

inf_int = data.frame("Flat_All" = integral(testrun_flat, 3000), 
                     "Single_1" = integral(single_list[[1]], 3000),
                     "Single_2" = integral(single_list[[2]], 3000),
                     "Single_3" = integral(single_list[[3]], 3000),
                     "Diff_1" = integral(diff_tax_list[[1]], 3000),
                     "Diff_2" = integral(diff_tax_list[[2]], 3000),
                     "Diff_3" = integral(diff_tax_list[[3]], 3000),
                     "Diff_4" = integral(diff_tax_list[[4]], 3000),
                     "Diff_5" = integral(diff_tax_list[[5]], 3000),
                     "Diff_6" = integral(diff_tax_list[[6]], 3000))

rownames(inf_int) <- c("Total Infections", "Average Resistance")


# Plot Baseline -----------------------------------------------------------

OG_melt <- melt(as.matrix(inf_int), measure.vars = colnames(inf_int))
rescale_data <- t(apply(inf_int, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
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

# Rescaled for Denominator ------------------------------------------------


diff_tax_list <- list()

for(i in 1:6) {
  diff_tax_list[[i]] <- data.frame(multi_int_fun(i, 365*3, parms, init, amr, agg_func)[[2]])
}

#Find The Resistance Lost 
#Baseline

base <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms)))
base_totinf <- integral(base, 3000)[1]
base_avgres_int <- integral(base, 3000)[2]

#Lost Resistance Function 
#Supply a dataframe with time and the effective tax over time for class 1, 2 and 3.
#I need to extract the 1-6 data with the 

#Now I need to do this for all 6 scenarios. 
parms_flat <- parms; parms_flat[grep("eff_tax", names(parms), value =T)]  <- 0.5
parms_single1 <- parms; parms_single1[grep("eff_tax1", names(parms), value =T)]  <- 0.5
parms_single2 <- parms; parms_single2[grep("eff_tax2", names(parms), value =T)]  <- 0.5
parms_single3 <- parms; parms_single3[grep("eff_tax3", names(parms), value =T)]  <- 0.5

parm_list <- list(parms_flat,parms_single1,parms_single2,parms_single3,parms_dual12,parms_dual23,parms_dual13)

over_parms <- append(parm_list, diff_tax_list)
list_usage <- lapply(over_parms, usage_fun)

# Relative Reduction Integral Comparison --------------------------------------------

reduc_usage_vec <- sapply(1:length(list_usage), function(x) sum(list_usage[[x]][,6]))
rel_AvgRes <- base_avgres_int - inf_int[2,]
rel_totinf <- inf_int[1,] - base_totinf

comp_totinf <- rel_totinf/reduc_usage_vec
comp_res <- rel_AvgRes/reduc_usage_vec

data.frame.rel <- rbind(comp_totinf, comp_res); rownames(data.frame.rel) <- c("Total Infections", "Average Resistance")

rescale_data_rel <- t(apply(data.frame.rel, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
OG_melt_rel <- melt(as.matrix(data.frame.rel), measure.vars = colnames(data.frame.rel))
rescale_melt_rel <- melt(as.matrix(rescale_data_rel), measure.vars = colnames(rescale_data_rel))
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
        strip.text = element_blank(), legend.position="none")
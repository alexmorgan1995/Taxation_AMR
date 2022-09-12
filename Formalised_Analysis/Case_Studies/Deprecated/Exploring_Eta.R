library("deSolve"); library("parallel");library("deSolve"); library("parallel"); library("ggpubr"); library("reshape2")
rm(list=ls())

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

# Shannon Function --------------------------------------------------------

shan_func <- function(data, parms, thresh) {
  
  data$aggR1 <- data$R1 + data$R12 + data$R13 + data$R123
  data$aggR2 <- data$R2 + data$R12 + data$R23 + data$R123
  data$aggR3 <- data$R3 + data$R13 + data$R23 + data$R123
  #Determine the X% Thresholds that you want to be under
  thresholds <- unlist(data[parms[["t_n"]]-1, 11:13]*thresh)
  
  under_thresh <- sweep(data[data[,1] > parms[["t_n"]],][,11:13], 2, thresholds)
  
  #Calculate the number of days you are under said threshold
  under_50 <- c(nrow(under_thresh[under_thresh$aggR1 < 0,]), 
                nrow(under_thresh[under_thresh$aggR2 < 0,]), 
                nrow(under_thresh[under_thresh$aggR3 < 0,]))
  
  #Find the Sum and make each value proportionate to one another 
  prop_vec <- under_50 / sum(under_50)
  prop_vec <- prop_vec[prop_vec != 0]
  return(-sum(sapply(1:length(prop_vec), function(x) prop_vec[x]*log(prop_vec[x]))))
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

usage_fun <- function(parms, data){
  tot_inf <- rowSums(data[3001:10001, 3:10])
  
  usage = data.frame("time" = seq(0,7000),
                     "PopUsage1" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax1_", x)]])*parms[["sigma1"]], 365*3)),
                                     rep((1-parms[["eff_tax1_6"]])*parms[["sigma1"]], 7001-365*3*6))*tot_inf,
                     "PopUsage2" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax2_", x)]])*parms[["sigma2"]], 365*3)),
                                     rep((1-parms[["eff_tax2_6"]])*parms[["sigma2"]], 7001-365*3*6))*tot_inf,
                     "PopUsage3" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax3_", x)]])*parms[["sigma3"]], 365*3)),
                                     rep((1-parms[["eff_tax3_6"]])*parms[["sigma3"]], 7001-365*3*6))*tot_inf)
  usage[usage < 0] <- 0
  usage$totusage = rowSums(usage[2:4])
  usage$reduc_use = (parms[["sigma1"]] + parms[["sigma2"]] + parms[["sigma3"]])*tot_inf - usage$totusage
  return(usage)
}

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = c(lambda = 1/365*(2), 
          beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
          r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
          r_t = 1/7, eta_wr = 0.3, eta_rw = 0.04, 
          eta_rr = 0.01, eta_rrr = 0.01,  
          c1 = 0.95, c2 = 0.925, c3 = 0.85,
          c12 = 0.85, c13 = 0.825, c23 = 0.75,
          c123 = 0.7,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# Run the Baseline Model --------------------------------------------------


run <- multi_int_fun(6, 365*3, parms,  init = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
                                                R12 = 0, R13 = 0, R23 = 0,
                                                R123 = 0), amr, agg_func)[[1]]

agg_run <- agg_func(run)
m_agg_run <- melt(agg_run, id.vars = "time", measure.vars = colnames(agg_run[4:6]))


ggplot(m_agg_run, aes(time, value, color = variable)) + geom_line()

# Testing Combinations of Eta ---------------------------------------------

eta_data <- expand.grid(c(0.01, 0.05, 0.1, 1), c(0.01, 0.05, 0.1, 1)); colnames(eta_data) <- c("eta_rw", "eta_wr")

eta_list <- list()

for(i in 1:nrow(eta_data)) {
  parms1 <- parms
  parms1[c("eta_rw", "eta_wr")] <- eta_data[i,]
  
  run <- multi_int_fun(6, 365*3, parms1,  init = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
                                                  R12 = 0, R13 = 0, R23 = 0,
                                                  R123 = 0), amr, agg_func)[[1]]
  
  agg_run <- agg_func(run)
  m_agg_run <- melt(agg_run, id.vars = "time", measure.vars = colnames(agg_run[4:6]))
  
  eta_list[[i]] <- ggplot(m_agg_run, aes(time, value, color = variable)) + geom_line() +
    labs(title = paste0("eta_rw = ",  eta_data[i,1], " | eta_wr = ",  eta_data[i,2]))
}

ggarrange(plotlist = eta_list, nrow = 4, ncol = 4, common.legend = T, legend = "bottom")


# Eta Heatmap -------------------------------------------------------------

eta_data_heat <- expand.grid(seq(0, 0.5, by = 0.025), seq(0, 0.5, by = 0.05)); colnames(eta_data_heat) <- c("eta_rw", "eta_wr")
eta_list <- data.frame(matrix(NA, nrow = nrow(eta_data_heat), ncol = 4))

for(i in 1:nrow(eta_data_heat)) {
  parms1 <- parms
  parms1[c("eta_rw", "eta_wr")] <- eta_data_heat[i,]
  
  run <- multi_int_fun(6, 365*3, parms1,  init = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
                                                   R12 = 0, R13 = 0, R23 = 0,
                                                   R123 = 0), amr, agg_func)[[1]]
  shan <- shan_func(run, parms1, 0.75)
  agg_run <- agg_func(run)
  agg_run_integral <- sum(rowMeans(agg_run[agg_run$time > 3000,4:6]))
  
  eta_list[i,] <- c(eta_data_heat[i,1],
                    eta_data_heat[i,2],
                    agg_run_integral,
                    shan)
  
  print(i/ nrow(eta_data_heat))
}


ggplot(eta_list, aes(X1, X2, fill= X3)) + labs(x = "Conversion from R to WT", y = "Conversion from WT to R", fill = "Average Resistance") +
  geom_tile() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + theme(legend.position = "bottom")

ggplot(eta_list, aes(X1, X2, fill= X4)) + labs(x = "Conversion from R to WT", y = "Conversion from WT to R", fill = "Shannon Index") +
  geom_tile() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + theme(legend.position = "bottom")

# Testing Scenarios -------------------------------------------------------


init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = c(lambda = 1/365*(2), 
          beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
          r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
          r_t = 1/7, eta_wr = 0.3, eta_rw = 0.5, 
          eta_rr = 0.01, eta_rrr = 0.01,  
          c1 = 0.95, c2 = 0.925, c3 = 0.85,
          c12 = 0.85, c13 = 0.825, c23 = 0.75,
          c123 = 0.7,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

m_data <- multi_int_fun(6, 365*3, parms,  init = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
                                                    R12 = 0, R13 = 0, R23 = 0,
                                                    R123 = 0), amr, agg_func)[[1]]

data_agg <- agg_func(m_data)

melt_data <- melt(data_agg, id.vars = "time", measure.vars = colnames(data_agg)[4:6])

ggplot(melt_data, aes(time, value, color = variable)) + geom_line()

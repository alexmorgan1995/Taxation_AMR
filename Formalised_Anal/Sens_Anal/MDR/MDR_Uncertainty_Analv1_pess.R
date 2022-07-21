library("deSolve"); library("parallel")
rm(list=ls())

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

# Baseline Parms ----------------------------------------------------------

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
          PED1 = 0.75, PED2 = 0.5, PED3 = 0.25, 
          t_n = 3000, time_between = Inf, rho = 0.1, base_tax = 0.5)

# The Function ------------------------------------------------------------

low_parm <- c(0, #lambda
              0, #beta
              1/50, #r_wt
              1/50, #r_r
              1/50, #r_rr
              1/50, #r_rrr
              1/50, #r_t
              0, #eta_wr
              0, #eta_rw
              0, #eta_rr
              0, #eta_rrr
              0.5, #c1
              0.5, #c2
              0.5, #c3
              0.5, #c12
              0.5, #c13
              0.5, #c23
              0.5, #c123
              0, #rho
              0) #baseline tax

high_parm <- c(1/300, #lambda
              10, #beta
              1/5, #r_wt
              1/5, #r_r
              1/5, #r_rr
              1/5, #r_rrr
              1/5, #r_t
              2, #eta_wr
              2, #eta_rw
              2, #eta_rr
              2, #eta_rrr
              1, #c1
              1, #c2
              1, #c3
              1, #c12
              1, #c13
              1, #c23
              1, #c123
              1, #rho
              1) #baseline tax

#Creating the Parm Dataframe

parm_data <- data.frame(t(replicate(100000, runif(20, low_parm, high_parm))))
colnames(parm_data) <- c("lambda", "beta", "r_wt", "r_r", "r_rr", "r_rrr","r_t",
                         "eta_wr", "eta_rw", "eta_rr", "eta_rrr",
                         "c1", "c2", "c3", "c12", "c13", "c23", "c123",  
                         "rho", "base_tax")

parm_data[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")][x,]), decreasing = F)))

parm_data[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")][x,]), decreasing = T)))

append_data <- data.frame(matrix(0, nrow = 100000, ncol = 18),
                          matrix(0.25, nrow = 100000, ncol = 3),
                          matrix(1, nrow = 100000, ncol = 3),
                          t_n = 3000,
                          time_between = Inf)

parm_data_comb <- data.frame(parm_data, append_data)
colnames(parm_data_comb)[21:38] <- names(parms)[22:39]
colnames(parm_data_comb)[39:41] <- names(parms)[3:5]
colnames(parm_data_comb)[42:44] <- names(parms)[40:42]

# Creating the Parallel Montonicity Function ------------------------------

mono_func <- function(n, parms_frame, init, amr_ode, usage_fun, multi_int_fun, low_parm, high_parm, agg_func) {
  parms_base = parms_frame[n,]
  #Run Baseline
  run_base <- remNA_func(data.frame(ode(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, hmax = 1)))
  run_base_agg <- agg_func(run_base)
  values <- tail(run_base, 1)
  
  if(values[4] == 0 & values[5] == 0 & values[6] == 0) {
    while(values[4] == 0 & values[5] == 0 & values[6] == 0) {
      parms_base[c(1:20)] <- runif(20,low_parm, high_parm)
      
      parms_base[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")] <- sort(as.numeric(parms_base[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")]), decreasing = F)
      parms_base[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")] <- 
        sort(as.numeric(parms_base[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")]), decreasing = T)
      
      run_base <- remNA_func(data.frame(ode(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, hmax = 1)))
      run_base_agg <- agg_func(run_base)
      values <- tail(run_base, 1)
    }
  }
  
  run <- run_base[run_base[,1] > parms[["t_n"]],]
  run_base_agg <- run_base_agg[run_base_agg[,1] > parms[["t_n"]],]
   
  base_tot_inf <- signif(sum(run[3:10]), 5)
  base_int_res <- signif(sum(rowMeans(run_base_agg[4:6]), 5))
  
  #Need to calculate a different baseline for each scenario for antibiotic usage 
  store_vec_res <- c()
  store_vec_inf <- c()
  
  for(i in 1:10){
    parms = parms_base
    if(i == 1) {
      parms[grep("eff_tax", names(parms), value =T)]  <- parms[["base_tax"]]
      out <- remNA_func(data.frame(ode(y = init, func = amr_ode, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
    if(i >= 2 & i <= 4) {
      parms[grep(paste0("eff_tax", i-1), names(parms), value =T)]  <- parms[["base_tax"]]
      out <- remNA_func(data.frame(ode(y = init, func = amr_ode, times = seq(0, 10000), parms = parms, hmax = 1)))
      
    }
    if(i >= 5 & i <= 10) {
      
      diff <- multi_int_fun(i-4, 365*3, parms, init, amr_ode, agg_func)
      out <- diff[[1]]
      parms <- diff[[2]]
    }
    
    data_temp <- out[out[,1] > parms[["t_n"]],]
    data_temp_agg <- agg_func(data_temp)
    
    out_vec <- signif(c(sum(data_temp[3:10]),
                        sum(rowMeans(data_temp_agg[4:6]))),5)
    
    reduc_usage_vec <- sum(usage_fun(parms)[,6])
    
    store_vec_inf[i] <- (out_vec[1] - base_tot_inf)/reduc_usage_vec
    store_vec_res[i] <- (base_int_res - out_vec[2])/reduc_usage_vec
  }
  
  output <- c(store_vec_inf, store_vec_res, parms_base[c(1:20)] )
  names(output) <- c("flat_inf", "single1_inf", "single2_inf", "single3_inf", "diff1_inf", "diff2_inf", "diff3_inf", "diff4_inf", "diff5_inf", "diff6_inf", 
                     "flat_res", "single1_res", "single2_res", "single3_res", "diff1_res", "diff2_res", "diff3_res", "diff4_res", "diff5_res", "diff6_res", 
                     names(parms_base[c(1:20)]))
  return(output)
}


# Run the Model ----------------------------------------------------------

start_time <- Sys.time()

test <- mclapply(1:5000, 
                 FUN = mono_func, 
                 parms_frame = parm_data_comb, 
                 init = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
                          R12 = 0, R13 = 0, R23 = 0,
                          R123 = 0), 
                 amr_ode = amr, 
                 usage_fun = usage_fun,
                 multi_int_fun = multi_int_fun,
                 low_parm = low_parm,
                 high_parm = high_parm,
                 agg_func = agg_func,
                 mc.cores = 10)

print(test)

comb_data <- data.frame(do.call(rbind, test))

saveRDS(parm_data_comb, "/cluster/home/amorgan/Sens_Anal_Output/parmFULL_MDRv1_pess.RDS")
saveRDS(comb_data, "/cluster/home/amorgan/Sens_Anal_Output/comb_dataFULLMDRv1_pess.RDS")

end_time <- Sys.time()
print(end_time - start_time)

library("deSolve"); library("sensitivity"); library("parallel");library("deSolve"); library("parallel"); library("ggpubr"); library("reshape2")

rm(list=ls())

# ODEs --------------------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    if((sigma1 + sigma2 + sigma3) > 1) {
      sigma1 <- sigma1 / (sigma1 + sigma2 + sigma3 + sigma_dummy)
      sigma2 <- sigma2 / (sigma1 + sigma2 + sigma3 + sigma_dummy)
      sigma3 <- sigma3 / (sigma1 + sigma2 + sigma3 + sigma_dummy)
    }
    
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

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = c(lambda = 1/365*(2), 
          beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25, sigma_dummy = 0.25,
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

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:6]),
                      sum(rowMeans(data_temp[4:6]))),5)
  
  return(out_vec)
}

# Functions for the Outcome Measures --------------------------------------

#The first Stage is to create a wrapper function for the eFAST analysis 
#The purpose of this is to have a function which can be used to take in a model input and to output the exact outcome measure you are looking to explore. 

ode_function <- function(x, init, outcome) {
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms = c(lambda = x[z,1], 
              beta = x[z,2], sigma1 = x[z,3], sigma2 = x[z,4], sigma3 = x[z,5], sigma_dummy = runif(1, 0 ,1),
              r_wt = x[z,6], r_r = x[z,7],  r_rr = x[z,8],  r_rrr = x[z,9], 
              r_t = x[z,10], eta_wr = x[z,11], eta_rw = x[z,12], 
              eta_rr = x[z,13], eta_rrr = x[z,14],  
              c1 = x[z,15], c2 = x[z,16], c3 = x[z,17],
              c12 = x[z,18], c13 = x[z,19], c23 = x[z,20],
              c123 = x[z,21],
              eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
              eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
              eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
              eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
              eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
              eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
              PED1 = x[z,22], PED2 = x[z,23], PED3 = x[z,24], 
              t_n = 3000, time_between = Inf, rho = x[z,25], base_tax = x[z,26])

    run <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms))
    
    #Cleaning the ODE outpout to ensure that there are no NAs or values below 
    
    additions <- length(grep("R",colnames(run)))
    run$X[is.na(run$X)] <- 1
    run[3:(3+additions)][is.na(run[3:(3+additions)]) | run[3:(3+additions)] < 0 ] <- 0
    
    #Progress Bar
    #print(paste0(round(z/nrow(x), digits  = 4)*100,"%"))
    
    agg_run <- agg_func(run)
    
    return_vec[z] <- c(sum(run[agg_run$time > 3000, 3:10]),
                       sum(rowMeans(agg_func(run)[agg_run$time > 3000,4:6])))[outcome]

  }
  return(return_vec)
}

# Run the FAST -----------------------------------------------------------

factors <- c("lambda", "beta", "sigma1", "sigma2", "sigma3", "r_wt", "r_r",  "r_rr",  "r_rrr", "r_t", 
             "eta_wr", "eta_rw", "eta_rr", "eta_rrr",  
             "c1", "c2", "c3", "c12", "c13", "c23", "c123",
             "PED1", "PED2", "PED3", 
             "rho", "base_tax")


init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

#Run the FBD Model
testfbd <- fast99(model = ode_function, factors = factors, n = 200, 
                  q.arg = list(list(min=0.0001, max=0.1), #lambda
                               list(min=0.0001, max=10), #beta
                               list(min=0, max=1), #sigma1
                               list(min=0, max=1), #sigma2
                               list(min=0, max=1), #sigma3
                               list(min=0.0001, max = 1/0.5), #r_wt
                               list(min=0.0001, max=1/0.5), #r_r
                               list(min=0.0001, max=1/0.5), #r_rr
                               list(min=0.0001, max=1/0.5), #r_rrr
                               list(min=0.0001, max=1/0.5), #r_t
                               list(min=0.0001, max=1), #eta_wr
                               list(min=0.0001, max=1), #eta_rw
                               list(min=0.0001, max=1), #eta_rr
                               list(min=0.0001, max=1), #eta_rrr
                               list(min=0.5, max=1), #c1
                               list(min=0.5, max=1), #c2
                               list(min=0.5, max=1), #c3
                               list(min=0.5, max=1), #c12
                               list(min=0.5, max=1), #c13
                               list(min=0.5, max=1), #c23
                               list(min=0.5, max=1), #c123
                               list(min=0, max=3), #PED1
                               list(min=0, max=3), #PED2
                               list(min=0, max=3), #PED3
                               list(min=0, max=1), #rho
                               list(min=0, max=1)), #basetax
                  init = init,
                  outcome = 1)

#Run the Resistance Model
testres <- fast99(model = ode_function, factors = factors, n = 200, 
                  q.arg = list(list(min=0.0001, max=0.1), #lambda
                               list(min=0.0001, max=10), #beta
                               list(min=0, max=1), #sigma1
                               list(min=0, max=1), #sigma2
                               list(min=0, max=1), #sigma3
                               list(min=0.0001, max = 1/0.5), #r_wt
                               list(min=0.0001, max=1/0.5), #r_r
                               list(min=0.0001, max=1/0.5), #r_rr
                               list(min=0.0001, max=1/0.5), #r_rrr
                               list(min=0.0001, max=1/0.5), #r_t
                               list(min=0.0001, max=1), #eta_wr
                               list(min=0.0001, max=1), #eta_rw
                               list(min=0.0001, max=1), #eta_rr
                               list(min=0.0001, max=1), #eta_rrr
                               list(min=0.5, max=1), #c1
                               list(min=0.5, max=1), #c2
                               list(min=0.5, max=1), #c3
                               list(min=0.5, max=1), #c12
                               list(min=0.5, max=1), #c13
                               list(min=0.5, max=1), #c23
                               list(min=0.5, max=1), #c123
                               list(min=0, max=3), #PED1
                               list(min=0, max=3), #PED2
                               list(min=0, max=3), #PED3
                               list(min=0, max=1), #rho
                               list(min=0, max=1)), #basetax
                  init = init,
                  outcome = 2)

# Plot the Output ---------------------------------------------------------

saveRDS(testfbd, "/cluster/home/amorgan/Sens_Anal_Output/Sens_Anal_FBD.RDS")
saveRDS(testres, "/cluster/home/amorgan/Sens_Anal_Output/Sens_Anal_Res.RDS")


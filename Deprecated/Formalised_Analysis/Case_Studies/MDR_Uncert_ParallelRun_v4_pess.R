library("deSolve"); library("parallel")
rm(list=ls())

# ODEs --------------------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    sigma_base1 <- sigma1; sigma_base2 <- sigma2; sigma_base3 <- sigma3
    
    if(t > t_n) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,1]*PED[1,1]) + (eff_tax[2,1]*PED[2,1]) + (eff_tax[3,1]*PED[3,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,1]*PED[1,2]) + (eff_tax[2,1]*PED[2,2]) + (eff_tax[3,1]*PED[3,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,1]*PED[1,3]) + (eff_tax[2,1]*PED[2,3]) + (eff_tax[3,1]*PED[3,3])))
    }
    
    if(t > (t_n + time_between)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,2]*PED[1,1]) + (eff_tax[2,2]*PED[2,1]) + (eff_tax[3,2]*PED[3,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,2]*PED[1,2]) + (eff_tax[2,2]*PED[2,2]) + (eff_tax[3,2]*PED[3,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,2]*PED[1,3]) + (eff_tax[2,2]*PED[2,3]) + (eff_tax[3,2]*PED[3,3])))
    }
    
    if(t > (t_n + time_between*2)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,3]*PED[1,1]) + (eff_tax[2,3]*PED[2,1]) + (eff_tax[3,3]*PED[3,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,3]*PED[1,2]) + (eff_tax[2,3]*PED[2,2]) + (eff_tax[3,3]*PED[3,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,3]*PED[1,3]) + (eff_tax[2,3]*PED[2,3]) + (eff_tax[3,3]*PED[3,3])))
    }
    
    if(t > (t_n + time_between*3)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,4]*PED[1,1]) + (eff_tax[2,4]*PED[2,1]) + (eff_tax[3,4]*PED[3,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,4]*PED[1,2]) + (eff_tax[2,4]*PED[2,2]) + (eff_tax[3,4]*PED[3,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,4]*PED[1,3]) + (eff_tax[2,4]*PED[2,3]) + (eff_tax[3,4]*PED[3,3])))
    }
    
    if(t > (t_n + time_between*4)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,5]*PED[1,1]) + (eff_tax[2,5]*PED[2,1]) + (eff_tax[3,5]*PED[3,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,5]*PED[1,2]) + (eff_tax[2,5]*PED[2,2]) + (eff_tax[3,5]*PED[3,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,5]*PED[1,3]) + (eff_tax[2,5]*PED[2,3]) + (eff_tax[3,5]*PED[3,3])))
    }
    
    if(t > (t_n + time_between*5)) {
      sigma1 <- sigma_base1*(1 + ((eff_tax[1,6]*PED[1,1]) + (eff_tax[2,6]*PED[2,1]) + (eff_tax[3,6]*PED[3,1])))
      sigma2 <- sigma_base2*(1 + ((eff_tax[1,6]*PED[1,2]) + (eff_tax[2,6]*PED[2,2]) + (eff_tax[3,6]*PED[3,2])))
      sigma3 <- sigma_base3*(1 + ((eff_tax[1,6]*PED[1,3]) + (eff_tax[2,6]*PED[2,3]) + (eff_tax[3,6]*PED[3,3])))
    }
    
    sigma1 <- ifelse(sigma1 > 0, sigma1, 0)
    sigma2 <- ifelse(sigma2 > 0, sigma2, 0)
    sigma3 <- ifelse(sigma3 > 0, sigma3, 0)
    
    if((sigma1 + sigma2 + sigma3) > 1) {
      sigma1 <- sigma1/(sigma1 + sigma2 + sigma3)
      sigma2 <- sigma2/(sigma1 + sigma2 + sigma3)
      sigma3 <- sigma3/(sigma1 + sigma2 + sigma3)
    }
    
    dX = lambda - lambda*X - beta*X*(Wt + R1*c1 + R2*c2 + R3*c3 + R12*c12 + R13*c13 + R23*c23 + R123*c123) +
      r_wt*Wt*(1 - (sigma1 + sigma2 + sigma3)) + r_r*R1*sigma1 + r_r*R2*sigma2 + r_r*R3*sigma3 +
      r_rr*R12*(sigma1 + sigma2) + r_rr*R13*(sigma1 + sigma3) + r_rr*R23*(sigma2 + sigma3) + 
      r_rrr*R123*(sigma1 + sigma2 + sigma3) + r_t*(1-rho)*(Wt*(sigma1 + sigma2 + sigma3) + R1*(sigma2 + sigma3) + 
                                                             R2*(sigma1 + sigma3) + R3*(sigma1 + sigma2) + 
                                                             R12*sigma3 + R13*sigma2 + R23*sigma1)
    
    dWt = - lambda*Wt + beta*X*Wt - r_wt*Wt*(1 - (sigma1 + sigma2 + sigma3)) - r_t*Wt*(1-rho)*(sigma1 + sigma2 + sigma3) +
      eta_rw*(R1 + R2 + R3 + R12 + R13 + R23 + R123)*(1 - (sigma1 + sigma2 + sigma3)) - 
      eta_wr*Wt*rho*(sigma1 + sigma2 + sigma3)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma2 + sigma3)*R1 - r_r*sigma1*R1 - eta_rr*R1*rho*sigma2 -
      eta_rr*R1*rho*sigma3 - eta_rw*R1*(1 - (sigma1 + sigma2 + sigma3)) + eta_wr*rho*Wt*sigma1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma1 + sigma3)*R2 - r_r*sigma2*R2 - eta_rr*R2*rho*sigma1 -
      eta_rr*R2*rho*sigma3 - eta_rw*R2*(1 - (sigma1 + sigma2 + sigma3)) + eta_wr*rho*Wt*sigma2
    
    dR3 = - lambda*R3 + beta*X*R3*c3 - r_t*(1-rho)*(sigma1 + sigma2)*R3 - r_r*sigma3*R3 - eta_rr*R3*rho*sigma1 -
      eta_rr*R3*rho*sigma2 - eta_rw*R3*(1 - (sigma1 + sigma2 + sigma3)) + eta_wr*rho*Wt*sigma3
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_t*(1-rho)*sigma3*R12 - r_rr*(sigma1 + sigma2)*R12 - eta_rrr*R12*rho*sigma3 -
      eta_rw*R12*(1 - (sigma1 + sigma2 + sigma3)) + eta_rr*R1*rho*sigma2 + eta_rr*R2*rho*sigma1 
    
    dR13 = - lambda*R13 + beta*X*R13*c13 - r_t*(1-rho)*sigma2*R13 - r_rr*(sigma1 + sigma3)*R13 - eta_rrr*R13*rho*sigma2 -
      eta_rw*R13*(1 - (sigma1 + sigma2 + sigma3)) + eta_rr*R1*rho*sigma3 + eta_rr*R3*rho*sigma1 
    
    dR23 = - lambda*R23 + beta*X*R23*c23 - r_t*(1-rho)*sigma1*R23 - r_rr*(sigma2 + sigma3)*R23 - eta_rrr*R23*rho*sigma1 -
      eta_rw*R23*(1 - (sigma1 + sigma2 + sigma3)) + eta_rr*R2*rho*sigma3 + eta_rr*R3*rho*sigma2 
    
    dR123 = - lambda*R123 + beta*X*R123*c123 - r_rrr*(sigma1 + sigma2 + sigma3)*R123 - eta_rw*R123*(1 - (sigma1 + sigma2 + sigma3)) + 
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
  
  if(timing[1,1] != max(dataframe[,1])) {
    for(i in 1:n_data){
      dataframe[seq(timing[[1]]+2,max(dataframe[,1])+1),i+1] <- timing[i,i+1]
    }
  }
  
  dataframe[dataframe < 1e-10] <- 0
  return(dataframe)
}

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, time_between, parms, init, func, agg_func){
  
  parms["time_between"] <- time_between
  
  #First Run
  run_1rd <- agg_func(remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))))
  
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  low_char_1rd <- names(values_1rd)[which.min(values_1rd)]
  high_char_1rd <- names(values_1rd)[which.max(values_1rd)]
  med_char_1rd <- names(values_1rd)[setdiff(1:3, c(which.min(values_1rd), which.max(values_1rd)))]
  
  parms[["eff_tax"]][as.numeric(substr(low_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[low_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][as.numeric(substr(high_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[high_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][as.numeric(substr(med_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[med_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][parms[["eff_tax"]] < 0] <- 0
  
  if(int_gen > 1) {
    
    for(i in 1:(int_gen-1)) {
      run <- agg_func(remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*i)), parms = parms))))
      values <- tail(run, 1)[4:6]
      
      if(values[1] == 0 & values[2] == 0 & values[3] == 0) {
        parms[["eff_tax"]][,i] <- 0
      } else {
        
        low_char <- names(values)[which.min(values)]
        high_char <- names(values)[which.max(values)]
        med_char <- names(values)[setdiff(1:3, c(which.min(values), which.max(values)))]
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i+1):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/values_1rd[med_char_1rd])))
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i+1):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/values_1rd[med_char_1rd])))
        parms[["eff_tax"]][as.numeric(substr(med_char, 2, 2)), c((i+1):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[med_char],1)/values_1rd[med_char_1rd])))
        parms[["eff_tax"]][parms[["eff_tax"]] < 0] <- 0
      }
    }
    
  }
  
  out_run <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
  
  return(list(out_run, parms))
}

# Single Tax Function -----------------------------------------------------

single_tax <- function(res_order, tax, parms, init, func, agg_func) {
  
  #First Run
  parms[["base_tax"]] <- tax
  
  run_1rd <- agg_func(remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))))
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  res_order_vec <- c(names(values_1rd)[which.max(values_1rd)],
                     names(values_1rd)[setdiff(1:3, c(which.min(values_1rd), which.max(values_1rd)))],
                     names(values_1rd)[which.min(values_1rd)])[res_order]
  
  parms[["eff_tax"]][as.numeric(substr(res_order_vec, 2, 2)), c(1:6)] <- as.numeric(parms[["base_tax"]])
  print(parms)
  #Real Model Run 
  run_real <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms)))
  return(run_real)
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
                     
                     "PopUsage1" = c(sapply(1:6, function(x) 
                       rep(parms[["sigma1"]] * (1 + ((parms[["eff_tax"]][1,x]*parms[["PED"]][1,1]) + (parms[["eff_tax"]][2,x]*parms[["PED"]][2,1]) + (parms[["eff_tax"]][3,x]*parms[["PED"]][3,1]))), 365*3)),
                       parms[["sigma1"]] * rep((1 + ((parms[["eff_tax"]][1,6]*parms[["PED"]][1,1]) + (parms[["eff_tax"]][2,6]*parms[["PED"]][2,1]) + (parms[["eff_tax"]][3,6]*parms[["PED"]][3,1]))), 
                                               7001-365*3*6)),
                     
                     "PopUsage2" = c(sapply(1:6, function(x) 
                       rep(parms[["sigma2"]] * (1 + ((parms[["eff_tax"]][1,x]*parms[["PED"]][1,2]) + (parms[["eff_tax"]][2,x]*parms[["PED"]][2,2]) + (parms[["eff_tax"]][3,x]*parms[["PED"]][3,2]))), 365*3)),
                       parms[["sigma2"]] * rep((1 + ((parms[["eff_tax"]][1,6]*parms[["PED"]][1,2]) + (parms[["eff_tax"]][2,6]*parms[["PED"]][2,2]) + (parms[["eff_tax"]][3,6]*parms[["PED"]][3,2]))), 
                                               7001-365*3*6)),
                     
                     "PopUsage3" = c(sapply(1:6, function(x) 
                       rep(parms[["sigma3"]] * (1 + ((parms[["eff_tax"]][1,x]*parms[["PED"]][1,3]) + (parms[["eff_tax"]][2,x]*parms[["PED"]][2,3]) + (parms[["eff_tax"]][3,x]*parms[["PED"]][3,3]))), 365*3)),
                       parms[["sigma3"]] * rep((1 + ((parms[["eff_tax"]][1,6]*parms[["PED"]][1,3]) + (parms[["eff_tax"]][2,6]*parms[["PED"]][2,3]) + (parms[["eff_tax"]][3,6]*parms[["PED"]][3,3]))), 
                                               7001-365*3*6)))
  usage[usage < 0] <- 0
  usage$totusage = rowSums(usage[2:4])
  usage$reduc_use = parms[["sigma1"]] + parms[["sigma2"]] + parms[["sigma3"]] - usage$totusage
  return(usage)
}

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = list(lambda = 1/365*(2), 
             beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
             r_t = 1/7, eta_wr = 0.3, eta_rw = 0.04, 
             eta_rr = 0.01, eta_rrr = 0.01,  
             c1 = 0.945, c2 = 0.925, c3 = 0.85,
             c12 = 0.845, c13 = 0.825, c23 = 0.75,
             c123 = 0.7,
             PED = matrix(c(-1.5, 0.4, 0.4, 
                            0.4, -1.25, 0.4,
                            0.4, 0.4, -1), #Be aware of this matrix
                          nrow = 3, ncol = 3, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 3, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# The Function ------------------------------------------------------------

low_parm <- c(1/3650*(2), #lambda
              0, #beta
              0, #sigma1
              0, #sigma2
              0, #sigma3
              1/50, #r_wt
              1/50, #r_r
              1/50, #r_rr
              1/50, #r_rrr
              1/50, #r_t
              0.03, #eta_wr
              0.004, #eta_rw
              0.001, #eta_rr
              0.001, #eta_rrr
              0.5, #c1
              0.5, #c2
              0.5, #c3
              0.5, #c12
              0.5, #c13
              0.5, #c23
              0.5, #c123
              0, #rho
              0) #baseline tax

high_parm <- c(1/36.5*(2), #lambda
               10, #beta
               1, #sigma1
               1, #sigma2
               1, #sigma3
               1/2, #r_wt
               1/2, #r_r
               1/2, #r_rr
               1/2, #r_rrr
               1/2, #r_t
               3, #eta_wr
               0.4, #eta_rw
               0.1, #eta_rr
               0.1, #eta_rrr
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

parm_data <- data.frame(t(replicate(10000, runif(23, low_parm, high_parm))))

colnames(parm_data) <- c("lambda", "beta", "sigma1", "sigma2", "sigma3", 
                         "r_wt", "r_r", "r_rr", "r_rrr","r_t",
                         "eta_wr", "eta_rw", "eta_rr", "eta_rrr",
                         "c1", "c2", "c3", "c12", "c13", "c23", "c123",  
                         "rho", "base_tax")

for(i in 1:nrow(parm_data)) {
  if(sum(parm_data[c("sigma1", "sigma2", "sigma3")][i,]) > 1) {
    parm_data[c("sigma1", "sigma2", "sigma3")][i,] <- parm_data[c("sigma1", "sigma2", "sigma3")][i,]/
      (sum(parm_data[c("sigma1", "sigma2", "sigma3")][i,]) + runif(1, 0, 1))
  }
}

parm_data[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr")][x,]), decreasing = T)))

parm_data[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")][x,]), decreasing = F)))


parm_data[c("c1", "c2", "c3")] <- t(sapply(1:nrow(parm_data), function(x) 
  sample(sort(as.numeric(parm_data[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")][x,]), decreasing = T)[1:3], 
         size = 3, replace = FALSE)))

parm_data[c("c12", "c13", "c23")] <-  t(sapply(1:nrow(parm_data), function(x) 
  sample(sort(as.numeric(parm_data[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")][x,]), decreasing = T)[4:6], 
         size = 3, replace = FALSE)))

parm_data["c123"] <- sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")][x,]), decreasing = T)[7])

parm_data_comb <- data.frame(parm_data, t_n = 3000,
                             time_between = Inf)

# Creating the Parallel Montonicity Function ------------------------------

mono_func <- function(n, parms_frame, init, amr_ode, usage_fun, multi_int_fun, low_parm, high_parm, agg_func, thresh) {
  
  parms_base = as.list(parms_frame[n,])
  parms_base = append(parms_base, parms["PED"]); parms_base = append(parms_base, parms["eff_tax"])
  
  #Run Baseline
  run_base <- remNA_func(data.frame(ode(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, hmax = 1)))
  run_base_agg <- agg_func(run_base)
  values <- tail(run_base_agg, 1)
  
  if(values[4] == 0 & values[5] == 0 & values[6] == 0) {
    while(values[4] == 0 & values[5] == 0 & values[6] == 0) {
      parms_base[c(1:23)] <- as.list(runif(23, low_parm, high_parm))
      
      if(sum(unlist(parms_base[c("sigma1", "sigma2", "sigma3")])) > 1) {
        parms_base[c("sigma1", "sigma2", "sigma3")] <- as.list(unlist(parms_base[c("sigma1", "sigma2", "sigma3")])/
                                                                 (sum(unlist(parms_base[c("sigma1", "sigma2", "sigma3")])) + runif(1, 0, 1)))
      }
      
      parms_base[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr")] <- as.list(sort(as.numeric(parms_base[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr")]), decreasing = T))
      parms_base[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")] <- as.list(sort(as.numeric(parms_base[c("r_wt", "r_r", "r_rr", "r_rrr", "r_t")]), decreasing = F))
      
      
      parms_base[c("c1", "c2", "c3")] <- 
        as.list(sample(sort(as.numeric(parms_base[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")]), decreasing = T)[1:3]), size = 3, replace = F)
      
      parms_base[c("c12", "c13", "c23")] <- 
        as.list(sample(sort(as.numeric(parms_base[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")]), decreasing = T)[4:6]), size = 3, replace = F)
      
      parms_base["c123"] <- 
        as.list(sort(as.numeric(parms_base[c("c1", "c2", "c3", "c12", "c13", "c23", "c123")]), decreasing = T)[7])
      
      run_base <- remNA_func(data.frame(ode(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, hmax = 1)))
      run_base_agg <- agg_func(run_base)
      values <- tail(run_base_agg, 1)
    }
  }
  
  run <- run_base[run_base[,1] > parms_base[["t_n"]],]
  run_base_agg <- run_base_agg[run_base_agg[,1] > parms_base[["t_n"]],]
  
  #Identifying the order of the resistances
  res_order_vec <- c(names(values[4:6])[which.max(values[4:6])],
                     names(values[4:6])[setdiff(1:3, c(which.min(values[4:6]), which.max(values[4:6])))],
                     names(values[4:6])[which.min(values[4:6])])
  
  #Storing info for the integrals 
  base_tot_inf <- signif(sum(run[3:10]), 5)
  base_int_res <- signif(sum(rowMeans(run_base_agg[4:6]), 5))
  
  #Need to calculate a different baseline for each scenario for antibiotic usage 
  store_vec_res <- c()
  store_vec_inf <- c()
  store_vec_shan <- c()
  
  for(i in 1:10){
    parms = parms_base
    if(i == 1) {
      parms[["eff_tax"]][,] <- parms[["base_tax"]]
      out <- remNA_func(data.frame(ode(y = init, func = amr_ode, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
    if(i >= 2 & i <= 4) {
      parms[["eff_tax"]][as.numeric(substr(res_order_vec[i-1], 2, 2)), c(1:6)] <- parms[["base_tax"]]
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
    
    #Aggregation
    out$aggR1 <- out$R1 + out$R12 + out$R13 + out$R123
    out$aggR2 <- out$R2 + out$R12 + out$R23 + out$R123
    out$aggR3 <- out$R3 + out$R13 + out$R23 + out$R123
    
    #Determine the X% Thresholds that you want to be under
    thresholds <- unlist(out[parms[["t_n"]]-1, 11:13]*thresh)
    
    under_thresh <- sweep(out[out[,1] > parms[["t_n"]],][,11:13], 2, thresholds)
    
    #Calculate the number of days you are under said threshold
    under_50 <- c(nrow(under_thresh[under_thresh$aggR1 < 0,]), 
                  nrow(under_thresh[under_thresh$aggR2 < 0,]), 
                  nrow(under_thresh[under_thresh$aggR3 < 0,]))
    
    #Find the Sum and make each value proportionate to one another 
    prop_vec <- under_50 / sum(under_50)
    prop_vec <- prop_vec[prop_vec != 0]
    
    #Store Computation Vectors 
    store_vec_inf[i] <- (out_vec[1] - base_tot_inf)/reduc_usage_vec
    store_vec_res[i] <- (base_int_res - out_vec[2])/reduc_usage_vec
    store_vec_shan[i] <- -sum(sapply(1:length(prop_vec), function(x) prop_vec[x]*log(prop_vec[x])))
  }
  
  output <- c(store_vec_inf, store_vec_res, store_vec_shan, parms_base[c(1:27)])
  names(output) <- c("flat_inf", "singleHR_inf", "singleMR_inf", "singleLR_inf", "diff1_inf", "diff2_inf", "diff3_inf", "diff4_inf", "diff5_inf", "diff6_inf", 
                     "flat_res", "singleHR_res", "singleMR_res", "singleLR_res", "diff1_res", "diff2_res", "diff3_res", "diff4_res", "diff5_res", "diff6_res",
                     "flat_shan", "singleHR_shan", "singleMR_shan", "singleLR_shan", "diff1_shan", "diff2_shan", "diff3_shan", "diff4_shan", "diff5_shan", "diff6_shan",
                     names(parms_base[c(1:27)]))
  
  return(output)
}

# Run the Model ----------------------------------------------------------

start_time <- Sys.time()

test <- mclapply(1:1000, 
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
                 thresh = 0.5,
                 mc.cores = 10) 

#Combine the Output into a "normal" looking dataframe
comb_data <- data.frame(do.call(rbind, test))
comb_data_new <- data.frame(matrix(NA, nrow = nrow(comb_data), ncol = 30))

for(i in 1:nrow(comb_data)) {
  comb_data_new[i,] <- unlist(comb_data[i,1:30])
}

colnames(comb_data_new) <- colnames(comb_data)[1:30]

#Update the Parameter Set 
parm_data_comb_new <- parm_data_comb
parm_data_comb_new[1:nrow(comb_data),1:23] <- comb_data[,31:53]

parm_list <- list()

for(i in 1:nrow(parm_data_comb_new)) {
  p_list <- as.list(unlist(parm_data_comb_new[i,]))
  p_list <- append(p_list, parms["eff_tax"])
  p_list <- append(p_list, parms["PED"])
  parm_list[[i]] <- p_list
}

#Save the output

saveRDS(parm_list, "/cluster/home/amorgan/Sens_Anal_Output/MDR_run_parms_v4_pess.RDS")
saveRDS(comb_data_new, "/cluster/home/amorgan/Sens_Anal_Output/MDR_run_v4_pess.RDS")

end_time <- Sys.time()
print(end_time - start_time)

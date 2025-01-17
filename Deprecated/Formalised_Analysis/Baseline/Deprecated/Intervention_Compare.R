library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

# Integral Function -------------------------------------------------------

integral <- function(data, t_n, thresh){
  
  #Aggregate the Data into Resistance Classes
  data$aggR1 <- data$R1 + data$R12 + data$R13 + data$R123
  data$aggR2 <- data$R2 + data$R12 + data$R23 + data$R123
  data$aggR3 <- data$R3 + data$R13 + data$R23 + data$R123
  
  #Subset the Dataframe so that only results after the intervention start
  data_temp <- data[data[,1] > t_n,]
  
  #Determine the X% Thresholds that you want to be under
  thresholds <- unlist(data[t_n-1, 11:13]* thresh)
  under_thresh <- sweep(data[data[,1] > t_n,][,11:13], 2, thresholds)

  #Calculate the number of days you are under said threshold
  under_50 <- c(nrow(under_thresh[under_thresh$aggR1 < 0,]), 
                nrow(under_thresh[under_thresh$aggR2 < 0,]), 
                nrow(under_thresh[under_thresh$aggR3 < 0,]))

  #Find the Sum and make each value proportionate to one another 
  prop_vec <- under_50 / sum(under_50)
  prop_vec <- prop_vec[prop_vec != 0]

  #Output the Optimisation Criteria 
  out_vec <- signif(c(sum(data_temp[3:10]),
                      sum(rowMeans(data_temp[11:13])),
                      -sum(sapply(1:length(prop_vec), function(x) prop_vec[x]*log(prop_vec[x])))), 5)
  
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

usage_fun <- function(parms, data){
  tot_inf <- rowSums(data[3001:10001,3:10])
  
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

# Intervention Scenarios --------------------------------------------------

#Flat Tax
parms1 <- parms; parms1[grep(paste0("eff_tax"), names(parms1), value = TRUE)] <- 0.5
testrun_flat <- list(remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))))

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
  dat <- multi_int_fun(i, 365*3, parms, init, amr, agg_func)[[1]]
  diff_tax_list[[i]] <- data.frame(dat)
}

# Obtain the Integrals ----------------------------------------------------

list_scen <- unlist(list(testrun_flat, single_list, diff_tax_list), recursive = FALSE)

opt_crit = data.frame("Flat_All" = integral(list_scen[[1]], 3000, 0.75), 
                     "Single_1" = integral(list_scen[[2]], 3000, 0.75),
                     "Single_2" = integral(list_scen[[3]], 3000, 0.75),
                     "Single_3" = integral(list_scen[[4]], 3000, 0.75),
                     "Diff_1" = integral(list_scen[[5]], 3000, 0.75),
                     "Diff_2" = integral(list_scen[[6]], 3000, 0.75),
                     "Diff_3" = integral(list_scen[[7]], 3000, 0.75),
                     "Diff_4" = integral(list_scen[[8]], 3000, 0.75),
                     "Diff_5" = integral(list_scen[[9]], 3000, 0.75),
                     "Diff_6" = integral(list_scen[[10]], 3000, 0.75))

rownames(opt_crit) <- c("Tot_Inf", "AvgRes", "Shan_Ind")

rescale_data <- data.frame(t(apply(opt_crit[1:2,], MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
rescale_data["Shan_Ind",] <- unlist(c("Flat_All" = NA, (opt_crit[3,2:10] - min(opt_crit[3,2:10]))/diff(range(opt_crit[3,2:10]))))

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
  diff_tax_list[[i]] <- data.frame(multi_int_fun(i, 365*3, parms, init, amr, agg_func)[[2]], 
                                   "scen" = paste0("diff", i))
}

#Find The Resistance Lost 
#Baseline

base <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms))
base_totinf <- integral(base, 3000, 0.75)[1]
base_avgres_int <- integral(base, 3000, 0.75)[2]

#Lost Resistance Function 
#Supply a dataframe with time and the effective tax over time for class 1, 2 and 3.
#I need to extract the 1-6 data with the 

#Change in Usage From Baseline - Only for differential taxation models 

#Now I need to do this for all 6 scenarios. 
parms_flat <- parms; parms_flat[grep("eff_tax", names(parms), value =T)]  <- 0.5
parms_single1 <- parms; parms_single1[grep("eff_tax1", names(parms), value =T)]  <- 0.5
parms_single2 <- parms; parms_single2[grep("eff_tax2", names(parms), value =T)]  <- 0.5
parms_single3 <- parms; parms_single3[grep("eff_tax3", names(parms), value =T)]  <- 0.5

parm_list <- list(parms_flat,parms_single1,parms_single2,parms_single3)

over_parms <- append(parm_list, diff_tax_list)

list_usage <- list()

for(i in 1:10) {
  list_usage[[i]] <- usage_fun(over_parms[[i]], list_scen[[i]])
}

# Relative Reduction Integral Comparison --------------------------------------------

reduc_usage_vec <- sapply(1:length(list_usage), function(x) sum(list_usage[[x]][,6]))
rel_AvgRes <- base_avgres_int - opt_crit[2,]
rel_totinf <- opt_crit[1,] - base_totinf

comp_totinf <- rel_totinf/reduc_usage_vec
comp_res <- rel_AvgRes/reduc_usage_vec

data.frame.rel <- rbind(comp_totinf, comp_res); rownames(data.frame.rel) <- c("Total Infections", "Average Resistance")
data.frame.rel["Shannon Index",] <- opt_crit["Shan_Ind",]

rescale_data_rel <- data.frame(t(apply(data.frame.rel, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
rescale_data_rel["Shannon Index",] <- rescale_data["Shan_Ind",]

OG_melt_rel <- melt(as.matrix(data.frame.rel), measure.vars = colnames(data.frame.rel))
OG_melt_rel$Var2 <- factor(rep(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)", "Single Tax (LR)",
                               "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)", "Diff Tax (4 Rd)",
                               "Diff Tax (5 Rd)","Diff Tax (6 Rd)"), each = 3))

OG_melt_rel$Var2  <- factor(OG_melt_rel$Var2 , levels = unique(OG_melt_rel$Var2))

rescale_melt_rel <- melt(as.matrix(rescale_data_rel), measure.vars = colnames(rescale_data_rel))

rescale_melt_rel$Var2 <- factor(rep(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)", "Single Tax (LR)",
                               "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)", "Diff Tax (4 Rd)",
                               "Diff Tax (5 Rd)","Diff Tax (6 Rd)"), each = 3))

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


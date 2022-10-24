library("deSolve"); library("reshape2"); library("sensitivity")


rm(list=ls())

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms, sigma_use1, sigma_use2, sigma_use3) {
  with(as.list(c(y, parms)), {
    
    #Specify the time-varying functions
    
    sigma_use1 <- sigma_func1(t)
    sigma_use2 <- sigma_func2(t)
    sigma_use3 <- sigma_func3(t)
    
    #ODES Below
    
    dX = lambda - lambda*X - beta*X*(Wt + R1*c1 + R2*c2 + R3*c3 + R12*c12 + R13*c13 + R23*c23 + R123*c123) +
      r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + r_r*R1*(1-(sigma_use2 + sigma_use3)) + r_r*R2*(1-(sigma_use1 + sigma_use3)) + r_r*R3*(1-(sigma_use1 + sigma_use2)) +
      r_rr*R12*(1-sigma_use3) + r_rr*R13*(1-sigma_use2) + r_rr*R23*(1-sigma_use1) + 
      r_rrr*R123 + 
      r_t*(1-rho)*(Wt*(sigma_use1 + sigma_use2 + sigma_use3) + R1*(sigma_use2 + sigma_use3) + 
                     R2*(sigma_use1 + sigma_use3) + R3*(sigma_use1 + sigma_use2) + 
                     R12*sigma_use3 + R13*sigma_use2 + R23*sigma_use1)
    
    dWt = - lambda*Wt + beta*X*Wt - r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) - r_t*Wt*(1-rho)*(sigma_use1 + sigma_use2 + sigma_use3) +
      eta_rw*(R1 + R2 + R3 + R12 + R13 + R23 + R123)*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) - 
      eta_wr*Wt*rho*(sigma_use1 + sigma_use2 + sigma_use3)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma_use2 + sigma_use3)*R1 - r_r*(1-(sigma_use2 + sigma_use3))*R1 - eta_rr*R1*rho*sigma_use2 -
      eta_rr*R1*rho*sigma_use3 - eta_rw*R1*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma_use1 + sigma_use3)*R2 - r_r*(1-(sigma_use1 + sigma_use3))*R2 - eta_rr*R2*rho*sigma_use1 -
      eta_rr*R2*rho*sigma_use3 - eta_rw*R2*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use2
    
    dR3 = - lambda*R3 + beta*X*R3*c3 - r_t*(1-rho)*(sigma_use1 + sigma_use2)*R3 - r_r*(1-(sigma_use1 + sigma_use2))*R3 - eta_rr*R3*rho*sigma_use1 -
      eta_rr*R3*rho*sigma_use2 - eta_rw*R3*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use3
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_t*(1-rho)*sigma_use3*R12 - r_rr*(1-sigma_use3)*R12 - eta_rrr*R12*rho*sigma_use3 -
      eta_rw*R12*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R1*rho*sigma_use2 + eta_rr*R2*rho*sigma_use1 
    
    dR13 = - lambda*R13 + beta*X*R13*c13 - r_t*(1-rho)*sigma_use2*R13 - r_rr*(1-sigma_use2)*R13 - eta_rrr*R13*rho*sigma_use2 -
      eta_rw*R13*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R1*rho*sigma_use3 + eta_rr*R3*rho*sigma_use1 
    
    dR23 = - lambda*R23 + beta*X*R23*c23 - r_t*(1-rho)*sigma_use1*R23 - r_rr*(1-sigma_use1)*R23 - eta_rrr*R23*rho*sigma_use1 -
      eta_rw*R23*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R2*rho*sigma_use3 + eta_rr*R3*rho*sigma_use2 
    
    dR123 = - lambda*R123 + beta*X*R123*c123 - r_rrr*R123 - eta_rw*R123*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + 
      eta_rrr*rho*(sigma_use3*R12 + sigma_use2*R13 + sigma_use1*R23)
    
    return(list(c(dX,dWt,
                  dR1,dR2,dR3,
                  dR12,dR13,dR23,
                  dR123)))
  })
}

# Extract Sigmas for the ApproxFun Function -------------------------------

approx_sigma <- function(sigma_mat){
  
  usage = data.frame("time" = seq(0,10000),
                     "PopUsage1" = c(rep(sigma_mat[1,1], 3000),
                                     rep(sigma_mat[1,2], 365*3), rep(sigma_mat[1,3], 365*3), rep(sigma_mat[1,4], 365*3),
                                     rep(sigma_mat[1,5], 365*3), rep(sigma_mat[1,6], 365*3), rep(sigma_mat[1,7], 10001 - (3000 + (365*3)*5))),
                     
                     "PopUsage2" = c(rep(sigma_mat[2,1], 3000),
                                     rep(sigma_mat[2,2], 365*3), rep(sigma_mat[2,3], 365*3), rep(sigma_mat[2,4], 365*3),
                                     rep(sigma_mat[2,5], 365*3), rep(sigma_mat[2,6], 365*3), rep(sigma_mat[2,7], 10001 - (3000 + (365*3)*5))),
                     
                     "PopUsage3" = c(rep(sigma_mat[3,1], 3000),
                                     rep(sigma_mat[3,2], 365*3), rep(sigma_mat[3,3], 365*3), rep(sigma_mat[3,4], 365*3),
                                     rep(sigma_mat[3,5], 365*3), rep(sigma_mat[3,6], 365*3), rep(sigma_mat[3,7], 10001 - (3000 + (365*3)*5))))
  return(usage)
}

# ODE Wrapper Function ----------------------------------------------------

ode_wrapper <- function(times, y, parms, func, approx_sigma) {
  
  sigma_mat = matrix(c(rep(parms[["sigma1"]], 7),
                       rep(parms[["sigma2"]], 7),
                       rep(parms[["sigma3"]], 7)), 
                     nrow = 3, ncol = 7, byrow = T)
  eff_tax <- parms[["eff_tax"]]; PED <- parms[["PED"]]
  
  if(parms[["int_round"]] > 0 ) {
    for(i in 1:parms[["int_round"]]) {
      stor_sigma <- sigma_mat[,i]
      
      sigma_mat[,(i+1):7] = c(stor_sigma[1]*(1 + ((eff_tax[1,i]*PED[1,1]) + (eff_tax[2,i]*PED[2,1]) + (eff_tax[3,i]*PED[3,1]))),
                              stor_sigma[2]*(1 + ((eff_tax[1,i]*PED[1,2]) + (eff_tax[2,i]*PED[2,2]) + (eff_tax[3,i]*PED[3,2]))),
                              stor_sigma[3]*(1 + ((eff_tax[1,i]*PED[1,3]) + (eff_tax[2,i]*PED[2,3]) + (eff_tax[3,i]*PED[3,3]))))
      
      
      sigma_mat[,(i+1):7][sigma_mat[,(i+1)] < 0.01] <- 0.01
      
      if(colSums(sigma_mat)[i+1] > 1) {
        sigma_mat[,(i+1):7] <- sigma_mat[,i+1]/(sum(sigma_mat[,i+1])+0.01)
      }
    }
  }
  
  parms[["sigma_mat"]] <- sigma_mat
  
  sigma_data <- approx_sigma(sigma_mat)
  
  sigma_func1 <<- approxfun(sigma_data[,c(1,2)], rule = 2)
  sigma_func2 <<- approxfun(sigma_data[,c(1,3)], rule = 2)
  sigma_func3 <<- approxfun(sigma_data[,c(1,4)], rule = 2)
  
  #Run the model 
  out <- data.frame(ode(y = init, func = func, times = times, parms = parms))
  
  n_data <- ncol(out)-1
  
  timing <- t(sapply(1:n_data, function(x)  out[max(which(!is.na(out[,x+1]))),]))
  
  if(timing[1,1] != tail(times,1)) {
    for(i in 1:n_data){
      out[seq(timing[[1]]+2,tail(times,1)+1),i+1] <- timing[i,i+1]
    }
  }
  out[out < 1e-10] <- 0
  
  return(list(out, parms))
}

# Function to Aggregate Resistance ----------------------------------------

agg_func <- function(data) {
  agg_data <- data.frame("time" = data$time,
                         "Susc" = data$X,
                         "WT" = data$Wt, 
                         "R1" = data$R1 + data$R12 + data$R13 + data$R123,
                         "R2" = data$R2 + data$R12 + data$R23 + data$R123,
                         "R3" = data$R3 + data$R13 + data$R23 + data$R123)
  return(agg_data)
}

# Single Taxation Function ------------------------------------------------

single_tax <- function(res_order, tax, parms, init, func, agg_func, ode_wrapper, approx_sigma) {
  
  #First Run
  parms[["base_tax"]] <- tax
  
  run_1rd <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms, approx_sigma)[[1]])
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  res_order_vec <- c(names(values_1rd)[which.max(values_1rd)],
                     names(values_1rd)[setdiff(1:3, c(which.min(values_1rd), which.max(values_1rd)))],
                     names(values_1rd)[which.min(values_1rd)])[res_order]
  
  parms[["eff_tax"]][as.numeric(substr(res_order_vec, 2, 2)), c(1:6)] <- as.numeric(parms[["base_tax"]])
  parms[["int_round"]] <- 1
  
  #Real Model Run 
  run_real <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms, approx_sigma)
  return(run_real)
}

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, time_between, parms, init, func, agg_func, ode_wrapper, approx_sigma){
  
  parms["time_between"] <- time_between
  
  #First Run
  run_1rd <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms, approx_sigma)[[1]])
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  low_char_1rd <- names(values_1rd)[which.min(values_1rd)]
  high_char_1rd <- names(values_1rd)[which.max(values_1rd)]
  med_char_1rd <- names(values_1rd)[setdiff(1:3, c(which.min(values_1rd), which.max(values_1rd)))]
  
  parms[["eff_tax"]][as.numeric(substr(low_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[low_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][as.numeric(substr(high_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[high_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][as.numeric(substr(med_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[med_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][parms[["eff_tax"]] < 0.00001] <- 0
  
  #First Round of Diff Taxation
  
  parms[["int_round"]] <- 1
  
  if(int_gen > 1) {
    #All Rounds Above 1
    for(i in 2:int_gen) {
      parms[["int_round"]] <- i-1
      run <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*(i-1))), parms = parms, approx_sigma)[[1]])
      values <- tail(run, 1)[4:6]
      
      if(values[1] == 0 & values[2] == 0 & values[3] == 0) {
        parms[["eff_tax"]][,i] <- 0
      } else {
        
        low_char <- names(values)[which.min(values)]
        high_char <- names(values)[which.max(values)]
        med_char <- names(values)[setdiff(1:3, c(which.min(values), which.max(values)))]
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/values_1rd[med_char_1rd])))
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/values_1rd[med_char_1rd])))
        parms[["eff_tax"]][as.numeric(substr(med_char, 2, 2)), c((i):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[med_char],1)/values_1rd[med_char_1rd])))
        parms[["eff_tax"]][parms[["eff_tax"]] < 0.00001] <- 0
      }
      parms[["int_round"]] <- i
    }
  }
  out_run <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms, approx_sigma)
  return(out_run)
}

# Integral ----------------------------------------------------------------

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
  prop_shan <- under_50 / sum(under_50)
  prop_shan <- prop_shan[prop_shan != 0]
  
  prop_avganti <- sum(under_50) / 7000
  
  #Output the Optimisation Criteria 
  out_vec <- signif(c(sum(data_temp[3:10]),
                      sum(rowMeans(data_temp[11:13])),
                      prop_avganti,
                      -sum(sapply(1:length(prop_shan), function(x) prop_shan[x]*log(prop_shan[x])))), 5)
  
  return(out_vec)
}


# Extract Usage for Integrals -----------------------------------------------------------

usage_fun <- function(parms){
  
  #In case there is no intervention
  if(parms[["int_round"]] == 0) {
    usage = data.frame("time" = seq(0,7000),
                       "PopUsage1" = rep(parms[["sigma_mat"]][1,2], 7001),
                       
                       "PopUsage2" = rep(parms[["sigma_mat"]][2,2], 7001),
                       
                       "PopUsage3" = rep(parms[["sigma_mat"]][3,2], 7001))
  }
  
  #Intervention
  if(parms[["int_round"]] > 0) {
    usage = data.frame("time" = seq(0,7000),
                       "PopUsage1" = c(rep(parms[["sigma_mat"]][1,2] , 365*3), rep(parms[["sigma_mat"]][1,3] , 365*3), rep(parms[["sigma_mat"]][1,4] , 365*3),
                                       rep(parms[["sigma_mat"]][1,5], 365*3), rep(parms[["sigma_mat"]][1,6] , 365*3), rep(parms[["sigma_mat"]][1,7] , (365*3)+(7001-(365*3)*6))),
                       
                       "PopUsage2" = c(rep(parms[["sigma_mat"]][2,2] , 365*3), rep(parms[["sigma_mat"]][2,3] , 365*3), rep(parms[["sigma_mat"]][2,4] , 365*3),
                                       rep(parms[["sigma_mat"]][2,5], 365*3), rep(parms[["sigma_mat"]][2,6] , 365*3), rep(parms[["sigma_mat"]][2,7] , (365*3)+(7001-(365*3)*6))),
                       
                       "PopUsage3" = c(rep(parms[["sigma_mat"]][3,2] , 365*3), rep(parms[["sigma_mat"]][3,3] , 365*3), rep(parms[["sigma_mat"]][3,4] , 365*3),
                                       rep(parms[["sigma_mat"]][3,5], 365*3), rep(parms[["sigma_mat"]][3,6] , 365*3), rep(parms[["sigma_mat"]][3,7] , (365*3)+(7001-(365*3)*6))))
  }
  
  usage$totusage = rowSums(usage[2:4])
  usage$reduc_use = parms[["sigma1"]] + parms[["sigma2"]] + parms[["sigma3"]] - usage$totusage
  
  return(usage)
}


# Base Parameters ---------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms_base = list(lambda = 1/365*(2), int_round = 1, 
             beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
             r_t = 1/7, eta_wr = 0.3, eta_rw = 0.04, 
             eta_rr = 0.01, eta_rrr = 0.01,  
             c1 = 0.945, c2 = 0.925, c3 = 0.85,
             c12 = 0.845, c13 = 0.825, c23 = 0.75,
             c123 = 0.7,
             PED = matrix(c(-1, 0.6, 0.6, 
                            0.6, -1, 0.6,
                            0.6, 0.6, -1), #Be aware of this matrix
                          nrow = 3, ncol = 3, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 3, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)


# Functions for the Outcome Measures --------------------------------------

#The first Stage is to create a wrapper function for the eFAST analysis 
#The purpose of this is to have a function which can be used to take in a model input and to output the exact outcome measure you are looking to explore. 

ode_function <- function(x, init, parms_base, outcome, intervention, integral, 
                         ode_wrapper, approx_sigma, multi_int_fun, usage_fun) {
  return_vec <- c()
  
  for(z in 1:nrow(x)) {
    
    #Create the Parameters
    parms = c(lambda = x[z,1], 
              beta = x[z,2], sigma1 = x[z,3], sigma2 = x[z,4], sigma3 = x[z,5],
              r_wt = x[z,6], r_r = x[z,7],  r_rr = x[z,8],  r_rrr = x[z,9], 
              r_t = x[z,10], eta_wr = x[z,11], eta_rw = x[z,12], 
              eta_rr = x[z,13], eta_rrr = x[z,14],  
              c1 = x[z,15], c2 = x[z,16], c3 = x[z,17],
              c12 = x[z,18], c13 = x[z,19], c23 = x[z,20],
              c123 = x[z,21], base_tax = x[z,23],rho = x[z,22], 
              t_n = 3000, time_between = Inf, int_round = 1)
    
    #need to rearrange the parameters here - especially sigma 
    
    if(parms["sigma1"] + parms["sigma2"] + parms["sigma3"] > 1) {
      parms[c("sigma1", "sigma2", "sigma3")] <- as.list(unlist(parms[c("sigma1", "sigma2", "sigma3")])/
                                                               (sum(unlist(parms[c("sigma1", "sigma2", "sigma3")])) + runif(1, 0, 1)))
      
    }
    
    saveRDS(parms, "/cluster/home/amorgan/Sensitivity/parms.RDS")
    
    parms <- as.list(parms)
    parms = append(parms, parms_base["PED"])
    parms = append(parms, parms_base["eff_tax"])
    
    #Run Baseline
    run_base <- ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms, approx_sigma)[[1]]
    run_base_agg <- agg_func(run_base)
    values <- tail(run_base_agg, 1)
    
    run <- run_base[run_base[,1] > parms[["t_n"]],]
    run_base_agg <- run_base_agg[run_base_agg[,1] > parms[["t_n"]],]
    tail_agg <<- tail(run_base_agg, 1)
    
    if(round(tail_agg[[2]], digits = 4) != 1) {
      #Identifying the order of the resistances
      res_order_vec <- c(names(values[4:6])[which.max(values[4:6])],
                         names(values[4:6])[setdiff(1:3, c(which.min(values[4:6]), which.max(values[4:6])))],
                         names(values[4:6])[which.min(values[4:6])])
      
      #Storing baseline info for the integrals 
      base_tot_inf <- signif(sum(run[3:10]), 5)
      base_int_res <- signif(sum(rowMeans(run_base_agg[4:6]), 5))
      
      #Single Taxation
      parms[["eff_tax"]][as.numeric(substr(res_order_vec[2], 2, 2)), c(1:6)] <- parms[["base_tax"]]
      parms[["int_round"]] <- 1
      out_run_single <- ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms, approx_sigma)
      out_single <- out_run_single[[1]]; parms_single <- out_run_single[[2]]
      
      #Differential Taxation 
      diff <- multi_int_fun(6, 365*3, parms, init, amr, agg_func, ode_wrapper, approx_sigma)
      out_diff <- diff[[1]]
      parms_diff <- diff[[2]]
      
      #Integrals
      integral_single <- integral(out_single, 3000, 0.5)
      integral_diff <- integral(out_diff, 3000, 0.5)
      
      reduc_usage_vec_single <- sum(usage_fun(parms_single)[,6])
      reduc_usage_vec_diff <- sum(usage_fun(parms_diff)[,6])
      
      #Store Computation Vectors 
      single_usage <- list((integral_single[1] - base_tot_inf)/reduc_usage_vec_single,
                           (base_int_res - integral_single[2])/reduc_usage_vec_single,
                           integral_single[3],
                           integral_single[4])[[outcome]]
      
      diff_usage <- list((integral_diff[1] - base_tot_inf)/reduc_usage_vec_diff,
                         (base_int_res - integral_diff[2])/reduc_usage_vec_diff,
                         integral_diff[3],
                         integral_diff[4])[[outcome]]
      return_vec[z] <- list(single_usage, diff_usage)[[intervention]]
    } else{
      return_vec[z] <- NA
    }
 
    print(paste0(round(z/nrow(x), digits  = 4)*100,"%"))
  }

  #Progress Bar
  return(return_vec)
}

# Test Parms --------------------------------------------------------------

parms <- data.frame(readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/parms.RDS"))

ode_function(x = parms, 
             init = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
                       R12 = 0, R13 = 0, R23 = 0,
                       R123 = 0), 
             parms_base = parms_base, 
             outcome = 1, 
             intervention = 1, 
             integral = integral, 
             ode_wrapper = ode_wrapper, 
             approx_sigma = approx_sigma, 
             multi_int_fun = multi_int_fun, 
             usage_fun = usage_fun)


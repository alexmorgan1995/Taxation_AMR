library("deSolve"); library("parallel")
rm(list=ls())

# ODEs --------------------------------------------------------------------


amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    #Set up the matrix for the Sigmas
    
    sigma_use1 <- sigma_func1(t)
    sigma_use2 <- sigma_func2(t)
    
    dX = lambda - lambda*X - beta*X*(Wt + R1*c1 + R2*c2 + R12*c12) +
      r_wt*Wt*(1 - (sigma_use1 + sigma_use2)) + r_r*R1*(1-sigma_use2) + r_r*R2*(1-sigma_use1) + 
      r_rr*R12 +
      + r_t*(1-rho)*(Wt*(sigma_use1 + sigma_use2) + R1*(sigma_use2) + R2*(sigma_use1))
    
    dWt = - lambda*Wt + beta*X*Wt - r_wt*Wt*(1 - (sigma_use1 + sigma_use2)) - r_t*Wt*(1-rho)*(sigma_use1 + sigma_use2) +
      eta_rw*(R1 + R2 + R12)*(1 - (sigma_use1 + sigma_use2)) - 
      eta_wr*Wt*rho*(sigma_use1 + sigma_use2)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma_use2)*R1 - r_r*(1-sigma_use2)*R1 - eta_rr*R1*rho*sigma_use2 - 
      eta_rw*R1*(1 - (sigma_use1 + sigma_use2)) + eta_wr*rho*Wt*sigma_use1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma_use1)*R2 - r_r*(1-sigma_use1)*R2 - eta_rr*R2*rho*sigma_use1 -
      eta_rw*R2*(1 - (sigma_use1 + sigma_use2)) + eta_wr*rho*Wt*sigma_use2
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_rr*R12 -
      eta_rw*R12*(1 - (sigma_use1 + sigma_use2)) + eta_rr*R1*rho*sigma_use2 + eta_rr*R2*rho*sigma_use1 
    
    return(list(c(dX,dWt,
                  dR1,dR2,
                  dR12)))
  })
}


# Approx Sigma Function ---------------------------------------------------

approx_sigma <- function(sigma_mat){
  
  usage = data.frame("time" = seq(0,10000),
                     "PopUsage1" = c(rep(sigma_mat[1,1], 3000),
                                     rep(sigma_mat[1,2], 365*3), rep(sigma_mat[1,3], 365*3), rep(sigma_mat[1,4], 365*3),
                                     rep(sigma_mat[1,5], 365*3), rep(sigma_mat[1,6], 365*3), rep(sigma_mat[1,7], 10001 - (3000 + (365*3)*5))),
                     
                     "PopUsage2" = c(rep(sigma_mat[2,1], 3000),
                                     rep(sigma_mat[2,2], 365*3), rep(sigma_mat[2,3], 365*3), rep(sigma_mat[2,4], 365*3),
                                     rep(sigma_mat[2,5], 365*3), rep(sigma_mat[2,6], 365*3), rep(sigma_mat[2,7], 10001 - (3000 + (365*3)*5))))
                     
  return(usage)
}

# ODE Wrapper -------------------------------------------------------------

ode_wrapper <- function(times, y, parms, func, approx_sigma) {
  
  sigma_mat = matrix(c(rep(parms[["sigma1"]], 7),
                       rep(parms[["sigma2"]], 7)), 
                     nrow = 2, ncol = 7, byrow = T)
  
  eff_tax <- parms[["eff_tax"]]; PED <- parms[["PED"]]
  
  if(parms[["int_round"]] > 0 ) {
    for(i in 1:parms[["int_round"]]) {
      stor_sigma <- sigma_mat[,i]
      
      sigma_mat[,(i+1):7] = c(stor_sigma[1]*(1 + ((eff_tax[1,i]*PED[1,1]) + (eff_tax[2,i]*PED[2,1]))),
                              stor_sigma[2]*(1 + ((eff_tax[1,i]*PED[1,2]) + (eff_tax[2,i]*PED[2,2]))))
      
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
                         "R1" = data$R1 + data$R12,
                         "R2" = data$R2 + data$R12)
  return(agg_data)
}

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, time_between, parms, init, func, agg_func, ode_wrapper, approx_sigma){
  
  parms["time_between"] <- time_between
  
  #First Run
  run_1rd <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms, approx_sigma)[[1]])
  values_1rd <- tail(run_1rd, 1)[4:5]
  
  low_char_1rd <- names(values_1rd)[which.min(values_1rd)]
  high_char_1rd <- names(values_1rd)[which.max(values_1rd)]
  med_char_val <- mean(as.numeric(values_1rd))
  
  parms[["eff_tax"]][as.numeric(substr(low_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[low_char_1rd]/med_char_val)))
  parms[["eff_tax"]][as.numeric(substr(high_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[high_char_1rd]/med_char_val)))
  
  #First Round of Diff Taxation
  
  parms[["int_round"]] <- 1
  
  if(int_gen > 1) {
    #All Rounds Above 1
    for(i in 2:int_gen) {
      parms[["int_round"]] <- i-1
      run <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*(i-1))), parms = parms, approx_sigma)[[1]])
      values <- tail(run, 1)[4:5]
      
      if(values[1] == 0 & values[2] == 0) {
        parms[["eff_tax"]][,i] <- 0
      } else {
        
        low_char <- names(values)[which.min(values)]
        high_char <- names(values)[which.max(values)]
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i):6)] <- ((1 + as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/med_char_val)))) -  
                                                                               (1 + parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)])) 
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i):6)] <- ((1 + as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/med_char_val)))) -  
                                                                                 (1 + parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)]))
        
      }
      parms[["int_round"]] <- i
    }
  }
  out_run <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms, approx_sigma)
  return(out_run)
}

# Single Taxation Function ------------------------------------------------

single_tax <- function(res_order, tax, parms, init, func, agg_func, ode_wrapper, approx_sigma) {
  
  #First Run
  parms[["base_tax"]] <- tax
  
  run_1rd <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms, approx_sigma)[[1]])
  values_1rd <- tail(run_1rd, 1)[4:5]
  
  res_order_vec <- c(names(values_1rd)[which.max(values_1rd)],
                     names(values_1rd)[which.min(values_1rd)])[res_order]
  
  parms[["eff_tax"]][as.numeric(substr(res_order_vec, 2, 2)), c(1:6)] <- as.numeric(parms[["base_tax"]])
  parms[["int_round"]] <- 1
  
  #Real Model Run 
  run_real <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms, approx_sigma)
  return(run_real)
}

# Extract Usage -----------------------------------------------------------

usage_fun <- function(parms){
  
  #In case there is no intervention
  if(parms[["int_round"]] == 0) {
    usage = data.frame("time" = seq(0,7000),
                       "PopUsage1" = rep(parms[["sigma_mat"]][1,2], 7001),
                       "PopUsage2" = rep(parms[["sigma_mat"]][2,2], 7001))
  }
  #Intervention
  if(parms[["int_round"]] > 0) {
    usage = data.frame("time" = seq(0,7000),
                       "PopUsage1" = c(rep(parms[["sigma_mat"]][1,2] , 365*3), rep(parms[["sigma_mat"]][1,3] , 365*3), rep(parms[["sigma_mat"]][1,4] , 365*3),
                                       rep(parms[["sigma_mat"]][1,5], 365*3), rep(parms[["sigma_mat"]][1,6] , 365*3), rep(parms[["sigma_mat"]][1,7] , (365*3)+(7001-(365*3)*6))),
                       "PopUsage2" = c(rep(parms[["sigma_mat"]][2,2] , 365*3), rep(parms[["sigma_mat"]][2,3] , 365*3), rep(parms[["sigma_mat"]][2,4] , 365*3),
                                       rep(parms[["sigma_mat"]][2,5], 365*3), rep(parms[["sigma_mat"]][2,6] , 365*3), rep(parms[["sigma_mat"]][2,7] , (365*3)+(7001-(365*3)*6))))
  }
  usage$totusage = rowSums(usage[2:3])
  usage$reduc_use = parms[["sigma1"]] + parms[["sigma2"]] - usage$totusage
  return(usage)
}



# Ban Function ------------------------------------------------------------

ban_wrapper <- function(times, init, parms, func, approx_sigma, ban) {
  
  sigma_mat = matrix(c(rep(parms[["sigma1"]], 7),
                       rep(parms[["sigma2"]], 7)), 
                     nrow = 2, ncol = 7, byrow = T)
  ban_vec <- c(0,0)
  PED <- parms[["PED"]]
  
  sigma_data <- approx_sigma(sigma_mat)
  
  sigma_func1 <<- approxfun(sigma_data[,c(1,2)], rule = 2)
  sigma_func2 <<- approxfun(sigma_data[,c(1,3)], rule = 2)
  
  #Run Baseline
  run_1rd <- agg_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms)))
  values_1rd <- tail(run_1rd, 1)[4:5]
  
  res_order_vec <- c(names(values_1rd)[which.max(values_1rd)],
                     names(values_1rd)[which.min(values_1rd)])[ban]
  
  #Run the Tax Model
  ban_vec <- c(0,0,0)
  PED <- parms[["PED"]]
  
  if(ban >= 1 & ban <= 2) {
    ban_vec[as.numeric(substr(res_order_vec, 2, 2))] <- 1
    sigma_mat[,2:7] = c(sigma_mat[1,1]*(1 + ((ban_vec[1]*PED[1,1]) + (ban_vec[2]*PED[2,1]))), 
                        sigma_mat[2,1]*(1 + ((ban_vec[1]*PED[1,2]) + (ban_vec[2]*PED[2,2]))))
    sigma_mat[as.numeric(substr(res_order_vec, 2, 2)),2:7] <- 0
  }
  
  if(ban == "all") {
    sigma_mat[,2:7] = c(0,0)
  }
  
  for(i in 1:7) {
    if(colSums(sigma_mat)[i] > 1) {
      
      sigma_mat[,i] <- sigma_mat[,i]/(sum(sigma_mat[,i])+0.01)
    }   
  }
  
  parms[["sigma_mat"]] <- sigma_mat; sigma_data <- approx_sigma(sigma_mat)
  sigma_func1 <<- approxfun(sigma_data[,c(1,2)], rule = 2)
  sigma_func2 <<- approxfun(sigma_data[,c(1,3)], rule = 2)
  
  #Run the model 
  parms[["int_round"]] <- 1
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


# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0,
          R12 = 0)

parms = list(lambda = 1/365*(2), 
             beta = 4.918742, sigma1 = 0.25, sigma2 = 0.25, 
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,
             r_t = 1/7, eta_wr = 1.53141359, eta_rw = 0.06203388, 
             eta_rr = 0.09420535, 
             c1 = 0.95636319, c2 = 0.90284600, 
             c12 = 0.62569857, 
             PED = matrix(c(-1.5, 0.5,
                            0, -1), #Be aware of this matrix
                          nrow = 2, ncol = 2, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 2, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# The Function ------------------------------------------------------------

low_parm <- c(1/3650*(2), #lambda
              0, #beta
              0, #sigma1
              0, #sigma2
              1/365, #r_wt
              1/365, #r_r
              1/365, #r_rr
              1/365, #r_t
              0.15, #eta_wr
              0.006, #eta_rw
              0.009, #eta_rr
              0.5, #c1
              0.5, #c2
              0.5, #c12
              0, #rho
              0) #baseline tax

high_parm <- c(1/36.5*(2), #lambda
               10, #beta
               1, #sigma1
               1, #sigma2
               1/1, #r_wt
               1/1, #r_r
               1/1, #r_rr
               1/1, #r_t
               10, #eta_wr
               0.6, #eta_rw
               0.9, #eta_rr
               1, #c1
               1, #c2
               1, #c12
               1, #rho
               1) #baseline tax

#Creating the Parm Dataframe

parm_data <- data.frame(t(replicate(10000, runif(16, low_parm, high_parm))))

colnames(parm_data) <- c("lambda", "beta", "sigma1", "sigma2", 
                         "r_wt", "r_r", "r_rr","r_t",
                         "eta_wr", "eta_rw", "eta_rr",
                         "c1", "c2", "c12",
                         "rho", "base_tax")

for(i in 1:nrow(parm_data)) {
  if(sum(parm_data[c("sigma1", "sigma2")][i,]) > 1) {
    parm_data[c("sigma1", "sigma2")][i,] <- parm_data[c("sigma1", "sigma2")][i,]/
      (sum(parm_data[c("sigma1", "sigma2")][i,]) + runif(1, 0, 1))
  }
}

parm_data[c("eta_wr", "eta_rw", "eta_rr")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("eta_wr", "eta_rw", "eta_rr")][x,]), decreasing = T)))

parm_data[c("r_wt", "r_r", "r_rr","r_t")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("r_wt", "r_r", "r_rr", "r_t")][x,]), decreasing = F)))


parm_data[c("c1", "c2")] <- t(sapply(1:nrow(parm_data), function(x) 
  sample(sort(as.numeric(parm_data[c("c1", "c2", "c12")][x,]), decreasing = T)[1:2], 
         size = 2, replace = FALSE)))

parm_data["c12"] <- sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("c1", "c2", "c12")][x,]), decreasing = T)[3])

parm_data_comb <- data.frame(parm_data, t_n = 3000,int_round = 0,
                             time_between = Inf)

# Creating the Parallel Montonicity Function ------------------------------

mono_func <- function(n, parms_frame, init, amr_ode, usage_fun, multi_int_fun, low_parm, high_parm, agg_func, thresh, ode_wrapper, approx_sigma) {
  
  parms_base = as.list(parms_frame[n,])
  parms_base = append(parms_base, parms["PED"])
  parms_base = append(parms_base, parms["eff_tax"])
  parms_base = append(parms_base, parms["sigma_mat"])
  
  #Run Baseline
  run_base <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, approx_sigma)[[1]]
  run_base_agg <- agg_func(run_base)
  values <- tail(run_base_agg, 1)
  
  if(values[4] == 0 & values[5] == 0) {
    while(values[4] == 0 & values[5] == 0) {
      parms_base[c(1:16)] <- as.list(runif(16, low_parm, high_parm))
      
      if(sum(unlist(parms_base[c("sigma1", "sigma2")])) > 1) {
        parms_base[c("sigma1", "sigma2")] <- as.list(unlist(parms_base[c("sigma1", "sigma2")])/
                                                       (sum(unlist(parms_base[c("sigma1", "sigma2")])) + runif(1, 0, 1)))
      }
      parms_base[c("eta_wr", "eta_rw", "eta_rr")] <- as.list(sort(as.numeric(parms_base[c("eta_wr", "eta_rw", "eta_rr")]), decreasing = T))
      parms_base[c("r_wt", "r_r", "r_rr", "r_t")] <- as.list(sort(as.numeric(parms_base[c("r_wt", "r_r", "r_rr", "r_t")]), decreasing = F))
      parms_base[c("c1", "c2")] <- 
        as.list(sample(sort(as.numeric(parms_base[c("c1", "c2", "c12")]), decreasing = T)[1:2]), size = 3, replace = F)
      parms_base["c123"] <- as.list(sort(as.numeric(parms_base[c("c1", "c2", "c12")]), decreasing = T)[3])
      
      run_base <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, approx_sigma)[[1]]
      run_base_agg <- agg_func(run_base)
      values <- tail(run_base_agg, 1)
    }
  }
  
  run <- run_base[run_base[,1] > parms_base[["t_n"]],]
  run_base_agg <- run_base_agg[run_base_agg[,1] > parms_base[["t_n"]],]
  
  #Identifying the order of the resistances
  res_order_vec <- c(names(values[4:5])[which.max(values[4:5])],
                     names(values[4:5])[which.min(values[4:5])])
  
  #Storing info for the integrals 
  base_tot_inf <- signif(sum(run[3:5]), 5)
  base_int_res <- signif(sum(rowMeans(run_base_agg[4:5])), 5)
  
  #Need to calculate a different baseline for each scenario for antibiotic usage 
  store_vec_res <- c()
  store_vec_inf <- c()
  store_vec_shan <- c()
  store_vec_avganti <- c()
  
  for(i in 1:11){
    parms = parms_base
    if(i == 1) {
      parms[["eff_tax"]][,] <- parms[["base_tax"]]
      parms[["int_round"]] <- 1
      out_run <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10000), parms = parms, approx_sigma)
      out <- out_run[[1]]
      parms <- out_run[[2]]
    }
    if(i >= 2 & i <= 3) {
      parms[["eff_tax"]][as.numeric(substr(res_order_vec[i-1], 2, 2)), c(1:6)] <- parms[["base_tax"]]
      parms[["int_round"]] <- 1
      out_run <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10000), parms = parms, approx_sigma)
      out <- out_run[[1]]
      parms <- out_run[[2]]
    }
    if(i >= 4 & i <= 9) {
      diff <- multi_int_fun(i-3, 365*3, parms, init, amr_ode, agg_func, ode_wrapper, approx_sigma)
      out <- diff[[1]]
      parms <- diff[[2]]
    }
    if(i >= 10 & i <= 11) {
      ban <- ban_wrapper(times = seq(0, 10000), init, parms, amr_ode, approx_sigma, ban = i-9)
      out <- ban[[1]]
      parms <- ban[[2]]
    }
    
    data_temp <- out[out[,1] > parms[["t_n"]],]
    data_temp_agg <- agg_func(data_temp)
    
    out_vec <- signif(c(sum(data_temp[3:5]),
                        sum(rowMeans(data_temp_agg[4:5]))),5)
    reduc_usage_vec <- sum(usage_fun(parms)[,5])
    
    #Aggregation
    out$aggR1 <- out$R1 + out$R12 
    out$aggR2 <- out$R2 + out$R12 
    
    #Determine the X% Thresholds that you want to be under
    thresholds <- unlist(out[parms[["t_n"]]-1,7:8]*thresh)
    under_thresh <- sweep(out[out[,1] > parms[["t_n"]],][,7:8], 2, thresholds)
    
    #Calculate the number of days you are under said threshold
    under_50 <- c(nrow(under_thresh[under_thresh$aggR1 < 0,]), 
                  nrow(under_thresh[under_thresh$aggR2 < 0,]))
    
    #Find the Sum and make each value proportionate to one another 
    prop_vec <- sum(under_50) / (10000 - parms[["t_n"]])
    
    prop_vec_shan <- under_50 / sum(under_50)
    prop_vec_shan <- prop_vec_shan[prop_vec_shan != 0]
    
    #Store Computation Vectors 
    if((base_int_res - out_vec[2]) < 0 & reduc_usage_vec < 0) {
      store_vec_res[i] <- -1000
    } else {
      store_vec_res[i] <- (base_int_res - out_vec[2])/reduc_usage_vec
    }
    store_vec_inf[i] <- (base_tot_inf - out_vec[1])/reduc_usage_vec
    store_vec_shan[i] <- -sum(sapply(1:length(prop_vec_shan), function(x) prop_vec_shan[x]*log(prop_vec_shan[x])))
    store_vec_avganti[i] <- prop_vec
  }
  
  output <- c(store_vec_inf, store_vec_res, store_vec_shan, store_vec_avganti, parms_base[c(1:22)])
  names(output) <- c("flat_inf", "singleHR_inf", "singleLR_inf", "diff1_inf", "diff2_inf", "diff3_inf", "diff4_inf", "diff5_inf", "diff6_inf", "banHR_inf", "banLR_inf", 
                     "flat_res", "singleHR_res", "singleLR_res", "diff1_res", "diff2_res", "diff3_res", "diff4_res", "diff5_res", "diff6_res", "banHR_res", "banLR_res", 
                     "flat_shan", "singleHR_shan", "singleLR_shan", "diff1_shan", "diff2_shan", "diff3_shan", "diff4_shan", "diff5_shan", "diff6_shan", "banHR_shan",  "banLR_shan", 
                     "flat_avganti", "singleHR_avganti", "singleLR_avganti", "diff1_avganti", "diff2_avganti", "diff3_avganti", "diff4_avganti", "diff5_avganti", "diff6_avganti", "banHR_avganti", "banLR_avganti", 
                     names(parms_base[c(1:22)]))
  return(output)
}

# Run the Model ----------------------------------------------------------

start_time <- Sys.time()

test <- mclapply(1:1000, 
                 FUN = mono_func, 
                 parms_frame = parm_data_comb, 
                 init = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, 
                          R12 = 0), 
                 amr_ode = amr, 
                 usage_fun = usage_fun,
                 multi_int_fun = multi_int_fun,
                 low_parm = low_parm,
                 high_parm = high_parm,
                 agg_func = agg_func,
                 ode_wrapper = ode_wrapper, 
                 approx_sigma = approx_sigma,
                 thresh = 0.5,
                 mc.cores = 10) 

#Combine the Output into a "normal" looking dataframe
comb_data <- data.frame(do.call(rbind, test))
comb_data_new <- data.frame(matrix(NA, nrow = nrow(comb_data), ncol = 44))

for(i in 1:nrow(comb_data)) {
  comb_data_new[i,] <- unlist(comb_data[i,1:44])
}

colnames(comb_data_new) <- colnames(comb_data)[1:44]

#Update the Parameter Set 
parm_data_comb_new <- parm_data_comb
parm_data_comb_new[1:nrow(comb_data),1:19] <- comb_data[,45:63]

parm_list <- list()

for(i in 1:nrow(parm_data_comb_new)) {
  p_list <- as.list(unlist(parm_data_comb_new[i,]))
  p_list <- append(p_list, parms["eff_tax"])
  p_list <- append(p_list, parms["PED"])
  parm_list[[i]] <- p_list
}

#Save the output

saveRDS(parm_list, "/cluster/home/amorgan/Sens_Anal_Output/MDR_run_parms_two.RDS")
saveRDS(comb_data_new, "/cluster/home/amorgan/Sens_Anal_Output/MDR_run_two.RDS")

end_time <- Sys.time()
print(end_time - start_time)

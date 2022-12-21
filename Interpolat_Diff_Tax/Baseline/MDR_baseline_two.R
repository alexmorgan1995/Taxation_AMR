library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
library("rbenchmark")

rm(list=ls())

# ODEs --------------------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    #Set up the matrix for the Sigmas
    
    sigma_use1 <- sigma_mat[1,1]; sigma_use2 <- sigma_mat[2,1]
    
    #For One Round
    if(t > t_n & int_round > 0) { 
      sigma_use1 <- sigma_mat[1,2]; sigma_use2 <- sigma_mat[2,2]
    }
    
    if(t > (t_n + time_between) & int_round > 1) { 
      sigma_use1 <- sigma_mat[1,3]; sigma_use2 <- sigma_mat[2,3]
    }
    
    if(t > (t_n + time_between*2) & int_round > 2) { 
      sigma_use1 <- sigma_mat[1,4]; sigma_use2 <- sigma_mat[2,4]
    }
    
    if(t > (t_n + time_between*3) & int_round > 3) {
      sigma_use1 <- sigma_mat[1,5]; sigma_use2 <- sigma_mat[2,5]
    }
    
    if(t > (t_n + time_between*4) & int_round > 4) { 
      sigma_use1 <- sigma_mat[1,6]; sigma_use2 <- sigma_mat[2,6]
    }
    
    if(t > (t_n + time_between*5) & int_round > 5) {
      sigma_use1 <- sigma_mat[1,7]; sigma_use2 <- sigma_mat[2,7]
    }
    
    dX = lambda - lambda*X - beta*X*(Wt + R1*c1 + R2*c2 + R12*c12) +
      r_wt*Wt*(1 - (sigma_use1 + sigma_use2)) + r_r*R1*sigma_use1 + r_r*R2*sigma_use2 + 
      r_rr*R12*(sigma_use1 + sigma_use2) +
      + r_t*(1-rho)*(Wt*(sigma_use1 + sigma_use2) + R1*(sigma_use2) + R2*(sigma_use1))
    
    dWt = - lambda*Wt + beta*X*Wt - r_wt*Wt*(1 - (sigma_use1 + sigma_use2)) - r_t*Wt*(1-rho)*(sigma_use1 + sigma_use2) +
      eta_rw*(R1 + R2 + R12)*(1 - (sigma_use1 + sigma_use2)) - 
      eta_wr*Wt*rho*(sigma_use1 + sigma_use2)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma_use2)*R1 - r_r*sigma_use1*R1 - eta_rr*R1*rho*sigma_use2 - 
      eta_rw*R1*(1 - (sigma_use1 + sigma_use2)) + eta_wr*rho*Wt*sigma_use1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma_use1)*R2 - r_r*sigma_use2*R2 - eta_rr*R2*rho*sigma_use1 -
      eta_rw*R2*(1 - (sigma_use1 + sigma_use2)) + eta_wr*rho*Wt*sigma_use2
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_rr*(sigma_use1 + sigma_use2)*R12 -
      eta_rw*R12*(1 - (sigma_use1 + sigma_use2)) + eta_rr*R1*rho*sigma_use2 + eta_rr*R2*rho*sigma_use1 
    
    return(list(c(dX,dWt,
                  dR1,dR2,
                  dR12)))
  })
}

# Approx Sigma Function ---------------------------------------------------

approx_sigma <- function(sigma_mat){
  
  usage = data.frame("time" = seq(0,10300),
                     "PopUsage1" = c(rep(sigma_mat[1,1], 3000),
                                     rep(sigma_mat[1,2], 365*3), rep(sigma_mat[1,3], 365*3), rep(sigma_mat[1,4], 365*3),
                                     rep(sigma_mat[1,5], 365*3), rep(sigma_mat[1,6], 365*3), rep(sigma_mat[1,7], 10301 - (3000 + (365*3)*5))),
                     
                     "PopUsage2" = c(rep(sigma_mat[2,1], 3000),
                                     rep(sigma_mat[2,2], 365*3), rep(sigma_mat[2,3], 365*3), rep(sigma_mat[2,4], 365*3),
                                     rep(sigma_mat[2,5], 365*3), rep(sigma_mat[2,6], 365*3), rep(sigma_mat[2,7], 10301 - (3000 + (365*3)*5))))
  
  return(usage)
}

# ODE Function Wrapper ----------------------------------------------------

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
  parms[["actual_tax"]][,1:6] <- parms[["eff_tax"]][,1]
  #First Round of Diff Taxation
  
  parms[["int_round"]] <- 1
  
  if(int_gen > 1) {
    #All Rounds Above 1
    for(i in 2:int_gen) {
      parms[["int_round"]] <- i-1
      run <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*(i-1))), parms = parms, approx_sigma)[[1]])
      values <- tail(run, 1)[4:5]
      
      if(values[1] == 0 & values[2] == 0) {
        parms[["eff_tax"]][,i:6] <- 0; parms[["actual_tax"]][,i:6] <- 0
      } else {
        
        low_char <- names(values)[which.min(values)]
        high_char <- names(values)[which.max(values)]
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i):6)] <- (((1 + as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/med_char_val)))) -  
                                                                                (1 + parms[["actual_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)])) / 
                                                                               (1 + parms[["actual_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)]))
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i):6)] <- (((1 + as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/med_char_val)))) -  
                                                                                 (1 + parms[["actual_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)])) / 
                                                                                (1 + parms[["actual_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)]))
        
        #Store the Actual Tax Rate
        new_tax <- c((1 + as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/med_char_val)))) - 1,
                     (1 + as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/med_char_val)))) - 1)
        names(new_tax) <- c(low_char, high_char)
        parms[["actual_tax"]][,i:6] <- c(new_tax[["R1"]], new_tax[["R2"]])
      }
      parms[["int_round"]] <- i
    }
  }
  out_run <- ode_wrapper(y = init, func = func, times = seq(0, 10300), parms = parms, approx_sigma)
  return(out_run)
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
  run_real <- ode_wrapper(y = init, func = func, times = seq(0, 10300), parms = parms, approx_sigma)
  return(run_real)
}

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0,
          R12 = 0)

parms = list(lambda = 1/365*(2), int_round = 1, 
             beta = 5, sigma1 = 0.25, sigma2 = 0.25, 
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,
             r_t = 1/7, eta_wr = 0.3, eta_rw = 0.04, 
             eta_rr = 0.01, 
             c1 = 0.945, c2 = 0.8, 
             c12 = 0.7, 
             PED = matrix(c(-1, 0.4,
                            0.4, -1), #Be aware of this matrix
                          nrow = 2, ncol = 2, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 2, ncol = 6, byrow = T),
             actual_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 2, ncol = 6, byrow = T),
             sigma_mat = matrix(c(0, 0, 0, 0, 0, 0, 0, 
                                  0, 0, 0, 0, 0, 0, 0), 
                                nrow = 2, ncol = 7, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# Intervention Scenarios --------------------------------------------------

#Flat Tax

parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5; parms1[["int_round"]] <- 1
testrun_flat <- list(ode_wrapper(y = init, func = amr, times = seq(0, 10300), parms = parms1, approx_sigma)[[1]])

#Single Tax 
single_list <- list()

for(i in 1:2) {
  parms1 <- parms
  single_list[[i]] <- single_tax(i, 0.5, parms1, init, amr, agg_func, ode_wrapper, approx_sigma)[[1]]
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  diff_tax_list[[i]] <- multi_int_fun(i, 365*3, parms, init, amr, agg_func, ode_wrapper, approx_sigma)[[1]]
}

#Bans
ban_list <- list()
for(i in 1:2) {
  ban_list[[i]] <- ban_wrapper(times = seq(0, 10300), init, parms, amr, approx_sigma, ban = i)[[1]]
}

# Plotting the Scenarios --------------------------------------------------

#Create a combined list of all the scenarios
list_scen <- unlist(list(testrun_flat, single_list, diff_tax_list), recursive = FALSE)

melt_data <- list()

#Melt each one
for(i in 1:length(list_scen)) {
  data_agg <- agg_func(list_scen[[i]]) 
  colnames(data_agg)[4:5] <- c("High Res (HR)",  "Low Res (LR)") 
  melt_data[[i]] <- data.frame(melt(data_agg, id.vars = "time", measure.vars = colnames(data_agg)[4:5]),
                               "scen" = c("flat", "single1", "single2",
                                          "diff1", "diff2", "diff3", "diff4", "diff5", "diff6")[i])
}

p_data <- list()

#Plotting Loop
for(i in 1:length(melt_data)) {
  data <- melt_data[[i]]
  p_data[[i]] <- ggplot(data, aes(time, value, color = variable)) + geom_line() + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "Time", y = "Prevalence", title = c("Flat Tax",
                                                 "Single Tax (HR)","Single Tax (LR)",
                                                 "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)",
                                                 "Diff Tax (4 Rd)", "Diff Tax (5 Rd)", "Diff Tax (6 Rd)")[i], color = "")
}

ggarrange(p_data[[1]], "", "",
          p_data[[2]], p_data[[3]], "",
          p_data[[4]], p_data[[5]], p_data[[6]],
          p_data[[7]], p_data[[8]], p_data[[9]],
          labels = c("A", "", "",
                     "B", "", "",
                     "C", "", "",
                     "", "", ""), hjust = -.1,  nrow = 4, ncol = 3, common.legend = T, legend = "bottom")

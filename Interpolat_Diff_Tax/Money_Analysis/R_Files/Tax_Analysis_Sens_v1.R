library("deSolve"); library("reshape2"); library("parallel")
rm(list=ls())

# ODEs --------------------------------------------------------------------

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
  
  usage = data.frame("time" = seq(0,10300),
                     "PopUsage1" = c(rep(sigma_mat[1,1], 3000),
                                     rep(sigma_mat[1,2], 365*3), rep(sigma_mat[1,3], 365*3), rep(sigma_mat[1,4], 365*3),
                                     rep(sigma_mat[1,5], 365*3), rep(sigma_mat[1,6], 365*3), rep(sigma_mat[1,7], 10301 - (3000 + (365*3)*5))),
                     
                     "PopUsage2" = c(rep(sigma_mat[2,1], 3000),
                                     rep(sigma_mat[2,2], 365*3), rep(sigma_mat[2,3], 365*3), rep(sigma_mat[2,4], 365*3),
                                     rep(sigma_mat[2,5], 365*3), rep(sigma_mat[2,6], 365*3), rep(sigma_mat[2,7], 10301 - (3000 + (365*3)*5))),
                     
                     "PopUsage3" = c(rep(sigma_mat[3,1], 3000),
                                     rep(sigma_mat[3,2], 365*3), rep(sigma_mat[3,3], 365*3), rep(sigma_mat[3,4], 365*3),
                                     rep(sigma_mat[3,5], 365*3), rep(sigma_mat[3,6], 365*3), rep(sigma_mat[3,7], 10301 - (3000 + (365*3)*5))))
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
  run_real <- ode_wrapper(y = init, func = func, times = seq(0, 10300), parms = parms, approx_sigma)
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
  
  parms[["actual_tax"]][,1:6] <- parms[["eff_tax"]][,1]
  
  #First Round of Diff Taxation
  
  parms[["int_round"]] <- 1
  
  if(int_gen > 1) {
    #All Rounds Above 1
    for(i in 2:int_gen) {
      parms[["int_round"]] <- i-1
      run <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*(i-1))), parms = parms, approx_sigma)[[1]])
      values <- tail(run, 1)[4:6]
      
      if(values[1] == 0 & values[2] == 0 & values[3] == 0) {
        parms[["eff_tax"]][,i:6] <- 0; parms[["actual_tax"]][,i:6] <- 0
      } else {
        
        low_char <- names(values)[which.min(values)]
        high_char <- names(values)[which.max(values)]
        med_char <- names(values)[setdiff(1:3, c(which.min(values), which.max(values)))]
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i):6)] <- (((1+as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/values_1rd[med_char_1rd]))))-
                                                                                (1+parms[["actual_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)]))/
                                                                               (1+parms[["actual_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)]))
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i):6)] <- (((1+as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/values_1rd[med_char_1rd]))))-
                                                                                 (1+parms[["actual_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)]))/
                                                                                (1+parms[["actual_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)]))
        parms[["eff_tax"]][as.numeric(substr(med_char, 2, 2)), c((i):6)] <- (((1+as.numeric((parms[["base_tax"]]*(tail(run[med_char],1)/values_1rd[med_char_1rd]))))-
                                                                                (1+parms[["actual_tax"]][as.numeric(substr(med_char, 2, 2)), (i-1)]))/
                                                                               (1+parms[["actual_tax"]][as.numeric(substr(med_char, 2, 2)), (i-1)]))
        
        #Store the Actual Tax Rate
        new_tax <- c((1+as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/values_1rd[med_char_1rd])))) - 1,
                     (1+as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/values_1rd[med_char_1rd])))) - 1,
                     (1+as.numeric((parms[["base_tax"]]*(tail(run[med_char],1)/values_1rd[med_char_1rd])))) - 1)
        names(new_tax) <- c(low_char, high_char, med_char)
        parms[["actual_tax"]][,i:6] <- c(new_tax[["R1"]], new_tax[["R2"]], new_tax[["R3"]])
        
      }
      parms[["int_round"]] <- i
    }
  }
  out_run <- ode_wrapper(y = init, func = func, times = seq(0, 10300), parms = parms, approx_sigma)
  return(out_run)
}

# Convert DiffTax to EffTax --------------------------------------------------

conv_diff_efftax <- function(efftax_mat, gen) {
  efftax_mat <- efftax_mat + 1
  eff_tax <- matrix(NA, nrow = 3, ncol = 6)
  eff_tax[,1:6] <- efftax_mat[,1]
  if(gen > 1) {
    for(i in 2:gen) {
      eff_tax[,i:6] <- efftax_mat[,i] * eff_tax[,i-1]
    }
  }
  return(eff_tax-1)
}

# Creating the Time + EffTax Matrix ---------------------------------------

year_tax <- function(mat){
  
  time_data <- data.frame("time" = seq(1, 20),
                          "Tax1" =  c(rep(mat[1,1], 3), rep(mat[1,2] , 3), rep(mat[1,3] , 3),
                                      rep(mat[1,4], 3), rep(mat[1,5] , 3), rep(mat[1,6] , 5)),
                          
                          "Tax2" =  c(rep(mat[2,1], 3), rep(mat[2,2] , 3), rep(mat[2,3] , 3),
                                      rep(mat[2,4], 3), rep(mat[2,5] , 3), rep(mat[2,6] , 5)),
                          
                          "Tax3" =  c(rep(mat[3,1], 3), rep(mat[3,2] , 3), rep(mat[3,3] , 3),
                                      rep(mat[3,4], 3), rep(mat[3,5] , 3), rep(mat[3,6] , 5)))
  return(time_data)
}

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = list(lambda = 1/365*(2), int_round = 1, 
             beta = 4.918742, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
             r_t = 1/7, 
             eta_wr = 1.53141359, eta_rw = 0.06203388, 
             eta_rr = 0.09420535, eta_rrr = 0.09420535,  
             c1 = 0.95636319, c2 = 0.90284600, c3 = 0.66383335,
             c12 = 0.62569857, c13 = 0.59669175, c23 = 0.56935615,
             c123 = 0.54109666,
             PED = matrix(c(-1.5, 1, 0.5, 
                            0.5, -1.25, 0.75,
                            0.25, 0.5, -1), #Be aware of this matrix
                          nrow = 3, ncol = 3, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 3, ncol = 6, byrow = T),
             actual_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 3, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# Tax Data and Sales ------------------------------------------------------

sales <- read.csv("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/US_China_Group_Price.csv")
t_sales <- setNames(data.frame(t(sales[,-1])), sales[,1])

#Create LMIC Groupings 

LMIC_tax_vector <- c("G1" = sum(t_sales[t_sales$Group_China == 1,]$Revenue_LMIC), 
                     "G2" = sum(t_sales[t_sales$Group_China == 2,]$Revenue_LMIC), 
                     "G3" = sum(t_sales[t_sales$Group_China == 3,]$Revenue_LMIC))

#Create HIC Groupings

HIC_tax_vector <- c("G1" = sum(t_sales[t_sales$Group_US == 1,]$Revenue_HIC), 
                     "G2" = sum(t_sales[t_sales$Group_US == 2,]$Revenue_HIC), 
                     "G3" = sum(t_sales[t_sales$Group_US == 3,]$Revenue_HIC))

comb_taxvector <- list(HIC_tax_vector, LMIC_tax_vector)

# Parameter Ranges --------------------------------------------------------

low_parm <- c(1/3650*(2), #lambda
              0, #beta
              0, #sigma1
              0, #sigma2
              0, #sigma3
              1/365, #r_wt
              1/365, #r_r
              1/365, #r_rr
              1/365, #r_rrr
              1/365, #r_t
              0.15, #eta_wr
              0.006, #eta_rw
              0.009, #eta_rr
              0.009, #eta_rrr
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
               1/1, #r_wt
               1/1, #r_r
               1/1, #r_rr
               1/1, #r_rrr
               1/1, #r_t
               10, #eta_wr
               0.6, #eta_rw
               0.9, #eta_rr
               0.9, #eta_rrr
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

parm_data_comb <- data.frame(parm_data, t_n = 3000, int_round = 0,
                             time_between = Inf)

# Run Sensitivity Analysis ------------------------------------------------

mono_func <- function(n, parms_frame, init, amr_ode, usage_fun, multi_int_fun, low_parm, high_parm, agg_func, ode_wrapper, approx_sigma) {
  
  parms_base = as.list(parms_frame[n,])
  parms_base = append(parms_base, parms["PED"])
  parms_base = append(parms_base, parms["eff_tax"])
  parms_base = append(parms_base, parms["actual_tax"])
  
  #Run Baseline
  run_base <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10300), parms = parms_base, approx_sigma)[[1]]
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
      
      run_base <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10300), parms = parms_base, approx_sigma)[[1]]
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
  
  #Need to calculate a different baseline for each scenario for antibiotic usage 
  tax_vec <- c()
  
  for(i in 1:10){
    parms = parms_base
    if(i == 1) {
      parms[["eff_tax"]][,] <- parms[["base_tax"]]
      parms[["int_round"]] <- 1
      out_run <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10300), parms = parms, approx_sigma)
      parms <- out_run[[2]]
      eff_tax <- parms$eff_tax
    }
    if(i >= 2 & i <= 4) {
      parms[["eff_tax"]][as.numeric(substr(res_order_vec[i-1], 2, 2)), c(1:6)] <- parms[["base_tax"]]
      parms[["int_round"]] <- 1
      out_run <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10300), parms = parms, approx_sigma)
      parms <- out_run[[2]]
      eff_tax <- parms$eff_tax
    }
    if(i >= 5 & i <= 10) {
      diff <- multi_int_fun(i-4, 365*3, parms, init, amr_ode, agg_func, ode_wrapper, approx_sigma)
      parms <- diff[[2]]
      eff_tax <- parms$actual_tax
    }

    for(z in 1:2){
      data_year_new <- year_tax(eff_tax)
      tax_LMIC_HIC <- comb_taxvector[[z]]
        
      data_year_new[,2:4] <- sweep(data_year_new[,2:4], 2, tax_LMIC_HIC, "*")
      data_year_new$total <- rowSums(data_year_new[,2:4])
      tax_vec[seq(1,20, by = 2)[i]+(z-1)] <- sum(data_year_new$total)
    }
  }
  
  output <- tax_vec
  names(output) <- c("flat_HIC", "flat_LMIC",
                     "singleHR_HIC", "singleHR_LMIC",
                     "singleMR_HIC", "singleMR_LMIC",
                     "singleLR_HIC", "singleLR_LMIC", 
                     "diff1_HIC", "diff1_LMIC", 
                     "diff2_HIC", "diff2_LMIC",
                     "diff3_HIC", "diff3_LMIC", 
                     "diff4_HIC", "diff4_LMIC", 
                     "diff5_HIC", "diff5_LMIC", 
                     "diff6_HIC", "diff6_LMIC")
  return(output)
}

# Run the Model -----------------------------------------------------------

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
                 ode_wrapper = ode_wrapper,
                 approx_sigma = approx_sigma,
                 mc.cores = 10) 
 
# save --------------------------------------------------------------------

#Combine the Output into a "normal" looking dataframe
comb_data <- data.frame(do.call(rbind, test))
comb_data_new <- data.frame(matrix(NA, nrow = nrow(comb_data), ncol = 20))

for(i in 1:nrow(comb_data)) {
  comb_data_new[i,] <- unlist(comb_data[i,1:20])
}

colnames(comb_data_new) <- colnames(comb_data)[1:20]

#Save the output
saveRDS(comb_data_new, "/cluster/home/amorgan/Sens_Tax/taxlist_v1.RDS")

end_time <- Sys.time()
print(end_time - start_time)


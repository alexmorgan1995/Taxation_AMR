library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
library("rbenchmark")

rm(list=ls())

# ODEs --------------------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    #Set up the matrix for the Sigmas
    
    sigma_use1 <- sigma_func1(t)
    sigma_use2 <- sigma_func2(t)
    sigma_use3 <- sigma_func3(t)
    sigma_use4 <- sigma_func4(t)
    
    
    dX = lambda - lambda*X - 
      beta*X*(Wt + R1*c1 + R2*c2 + R3*c3 + R4*c4 + 
                R12*c12 + R13*c13 + R14*c14 + 
                R23*c23 + R24*c24 + R34*c34 + 
                R123*c123 + R124*c124 + R134*c134 + R234*c234 +
                R1234*c1234) +
      r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + 
      r_r*R1*(1-(sigma_use2 + sigma_use3 + sigma_use4)) + 
      r_r*R2*(1-(sigma_use1 + sigma_use3 + sigma_use4)) + 
      r_r*R3*(1-(sigma_use1 + sigma_use2 + sigma_use4)) + 
      r_r*R4*(1-(sigma_use1 + sigma_use2 + sigma_use3)) +
      r_rr*R12*(1-(sigma_use3 + sigma_use4)) + 
      r_rr*R13*(1-(sigma_use2 + sigma_use4)) + 
      r_rr*R23*(1-(sigma_use1 + sigma_use4)) + 
      r_rr*R14*(1-(sigma_use2 + sigma_use3)) + 
      r_rr*R24*(1-(sigma_use1 + sigma_use3)) + 
      r_rr*R34*(1-(sigma_use1 + sigma_use2)) + 
      r_rrr*R123*(1-(sigma_use4)) + 
      r_rrr*R124*(1-(sigma_use3)) + 
      r_rrr*R234*(1-(sigma_use1)) + 
      r_rrr*R134*(1-(sigma_use2)) + 
      r_rrrr*R1234 + 
      r_t*(1-rho)*(Wt*(sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4) + 
                     R1*(sigma_use2 + sigma_use3 + sigma_use4) + 
                     R2*(sigma_use1 + sigma_use3 + sigma_use4) + 
                     R3*(sigma_use1 + sigma_use2 + sigma_use4) + 
                     R4*(sigma_use1 + sigma_use2 + sigma_use3) + 
                     R12*(sigma_use3 + sigma_use4) + R13*(sigma_use2 + sigma_use4) + 
                     R14*(sigma_use2 + sigma_use3) + R23*(sigma_use1 + sigma_use4) +
                     R24*(sigma_use1 + sigma_use3) + R34*(sigma_use1 + sigma_use2) +
                     R123*(sigma_use4) + R124*(sigma_use3) + 
                     R234*(sigma_use1) + R134*(sigma_use2))
    
    dWt = - lambda*Wt + beta*X*Wt - 
      r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) - 
      r_t*Wt*(1-rho)*(sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4) +
      eta_rw*(R1 + R2 + R3 + R4 + 
                R12 + R13 + R14 + R23 + R24 + R34 +
                R123 + R124 + R134 + R234 + 
                R1234)*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) - 
      eta_wr*Wt*rho*(sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma_use2 + sigma_use3 + sigma_use4)*R1 - r_r*(1-(sigma_use2 + sigma_use3 + sigma_use4))*R1 - eta_rr*R1*rho*sigma_use2 -
      eta_rr*R1*rho*sigma_use3 - eta_rr*R1*rho*sigma_use4 - eta_rw*R1*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_wr*rho*Wt*sigma_use1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma_use1 + sigma_use3 + sigma_use4)*R2 - r_r*(1-(sigma_use1 + sigma_use3 + sigma_use4))*R2 - eta_rr*R2*rho*sigma_use1 -
      eta_rr*R2*rho*sigma_use3 - eta_rr*R2*rho*sigma_use4 - eta_rw*R2*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_wr*rho*Wt*sigma_use2
    
    dR3 = - lambda*R3 + beta*X*R3*c3 - r_t*(1-rho)*(sigma_use1 + sigma_use2 + sigma_use4)*R3 - r_r*(1-(sigma_use1 + sigma_use2 + sigma_use4))*R3 - eta_rr*R3*rho*sigma_use1 -
      eta_rr*R3*rho*sigma_use2 - eta_rr*R3*rho*sigma_use4 - eta_rw*R3*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_wr*rho*Wt*sigma_use3
    
    dR4 = - lambda*R4 + beta*X*R4*c4 - r_t*(1-rho)*(sigma_use1 + sigma_use2 + sigma_use3)*R4 - r_r*(1-(sigma_use1 + sigma_use2 + sigma_use3))*R4 - eta_rr*R4*rho*(sigma_use1 + sigma_use2 + sigma_use3) - 
      eta_rw*R4*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_wr*rho*Wt*sigma_use4
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_t*(1-rho)*(sigma_use3 + sigma_use4)*R12 - r_rr*(1-(sigma_use3 + sigma_use4))*R12 - eta_rrr*R12*rho*(sigma_use3) - eta_rrr*R12*rho*sigma_use4 -
      eta_rw*R12*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rr*R1*rho*sigma_use2 + eta_rr*R2*rho*sigma_use1 
    
    dR13 = - lambda*R13 + beta*X*R13*c13 - r_t*(1-rho)*(sigma_use2 + sigma_use4)*R13 - r_rr*(1-(sigma_use2 + sigma_use4))*R13 - eta_rrr*R13*rho*(sigma_use2) - eta_rrr*R13*rho*sigma_use4 -
      eta_rw*R13*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rr*R1*rho*sigma_use3 + eta_rr*R3*rho*sigma_use1 
    
    dR14 = - lambda*R14 + beta*X*R14*c14 - r_t*(1-rho)*(sigma_use2 + sigma_use3)*R14 - r_rr*(1-(sigma_use2 + sigma_use3))*R14 - eta_rrr*R14*rho*(sigma_use2) - eta_rrr*R14*rho*sigma_use3 -
      eta_rw*R14*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rr*R1*rho*sigma_use4 + eta_rr*R4*rho*sigma_use1 
    
    dR23 = - lambda*R23 + beta*X*R23*c23 - r_t*(1-rho)*(sigma_use1 + sigma_use4)*R23 - r_rr*(1-(sigma_use1 + sigma_use4))*R23 - eta_rrr*R23*rho*(sigma_use1) - eta_rrr*R23*rho*sigma_use4 -
      eta_rw*R23*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rr*R2*rho*sigma_use3 + eta_rr*R3*rho*sigma_use2
    
    dR24 = - lambda*R24 + beta*X*R24*c24 - r_t*(1-rho)*(sigma_use1 + sigma_use3)*R24 - r_rr*(1-(sigma_use1 + sigma_use3))*R24 - eta_rrr*R24*rho*(sigma_use1) - eta_rrr*R24*rho*sigma_use3 -
      eta_rw*R24*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rr*R2*rho*sigma_use4 + eta_rr*R4*rho*sigma_use2 
    
    dR34 = - lambda*R34 + beta*X*R34*c34 - r_t*(1-rho)*(sigma_use1 + sigma_use2)*R34 - r_rr*(1-(sigma_use1 + sigma_use2))*R34 - eta_rrr*R34*rho*(sigma_use1) - eta_rrr*R34*rho*sigma_use2 -
      eta_rw*R34*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rr*R3*rho*sigma_use4 + eta_rr*R4*rho*sigma_use3 
    
    
    dR123 = - lambda*R123 + beta*X*R123*c123 - r_t*(1-rho)*R123*(sigma_use4) - r_rrr*(1-(sigma_use4))*R123 - 
      eta_rw*R123*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rrr*rho*(sigma_use3*R12) + eta_rrr*rho*sigma_use2*R13 + eta_rrr*rho*sigma_use1*R23 - 
      eta_rrrr*rho*R123*(sigma_use4)
    
    dR124 = - lambda*R124 + beta*X*R124*c124 - r_t*(1-rho)*R124*(sigma_use3) - r_rrr*(1-(sigma_use3))*R124 - 
      eta_rw*R124*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rrr*rho*(sigma_use4*R12) + eta_rrr*rho*sigma_use2*R14 + eta_rrr*rho*sigma_use1*R24 - 
      eta_rrrr*rho*R124*(sigma_use3)
    
    dR134 = - lambda*R134 + beta*X*R134*c134 - r_t*(1-rho)*R134*(sigma_use2) - r_rrr*(1-(sigma_use2))*R134 - 
      eta_rw*R134*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rrr*rho*(sigma_use4*R13) + eta_rrr*rho*sigma_use1*R34 + eta_rrr*rho*sigma_use3*R14 - 
      eta_rrrr*rho*R134*(sigma_use2)
    
    dR234 = - lambda*R234 + beta*X*R234*c234 - r_t*(1-rho)*R234*(sigma_use1) - r_rrr*(1-(sigma_use1))*R234 - 
      eta_rw*R234*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + eta_rrr*rho*(sigma_use4*R23) + eta_rrr*rho*sigma_use2*R34 + eta_rrr*rho*sigma_use3*R24 - 
      eta_rrrr*rho*R234*(sigma_use1)
    
    dR1234 = - lambda*R1234 + beta*X*R1234*c1234 - r_rrrr*R1234 - eta_rw*R1234*(1 - (sigma_use1 + sigma_use2 + sigma_use3 + sigma_use4)) + 
      eta_rrrr*rho*(R123*(sigma_use4) + R124*(sigma_use3) + R134*(sigma_use2) + R234*(sigma_use1))
    
    return(list(c(dX,dWt,
                  dR1,dR2,dR3,dR4,
                  dR12,dR13,dR14,
                  dR23,dR24,dR34,
                  dR123,dR124,dR134,dR234,
                  dR1234)))
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
                                     rep(sigma_mat[3,5], 365*3), rep(sigma_mat[3,6], 365*3), rep(sigma_mat[3,7], 10001 - (3000 + (365*3)*5))),
                     
                     "PopUsage3" = c(rep(sigma_mat[4,1], 3000),
                                     rep(sigma_mat[4,2], 365*3), rep(sigma_mat[4,3], 365*3), rep(sigma_mat[4,4], 365*3),
                                     rep(sigma_mat[4,5], 365*3), rep(sigma_mat[4,6], 365*3), rep(sigma_mat[4,7], 10001 - (3000 + (365*3)*5))))
  return(usage)
}

# ODE Function Wrapper ----------------------------------------------------

ode_wrapper <- function(times, y, parms, func, approx_sigma) {
  
  sigma_mat = matrix(c(rep(parms[["sigma1"]], 7),
                       rep(parms[["sigma2"]], 7),
                       rep(parms[["sigma3"]], 7),
                       rep(parms[["sigma4"]], 7)), 
                     nrow = 4, ncol = 7, byrow = T)
  
  eff_tax <- parms[["eff_tax"]]; PED <- parms[["PED"]]
  
  if(parms[["int_round"]] > 0 ) {
    for(i in 1:parms[["int_round"]]) {
      stor_sigma <- sigma_mat[,i]
      
      sigma_mat[,(i+1):7] = c(stor_sigma[1]*(1 + ((eff_tax[1,i]*PED[1,1]) + (eff_tax[2,i]*PED[2,1]) + (eff_tax[3,i]*PED[3,1]) + (eff_tax[4,i]*PED[4,1]))),
                              stor_sigma[2]*(1 + ((eff_tax[1,i]*PED[1,2]) + (eff_tax[2,i]*PED[2,2]) + (eff_tax[3,i]*PED[3,2]) + (eff_tax[4,i]*PED[4,2]))),
                              stor_sigma[3]*(1 + ((eff_tax[1,i]*PED[1,3]) + (eff_tax[2,i]*PED[2,3]) + (eff_tax[3,i]*PED[3,3]) + (eff_tax[4,i]*PED[4,3]))),
                              stor_sigma[4]*(1 + ((eff_tax[1,i]*PED[1,4]) + (eff_tax[2,i]*PED[2,4]) + (eff_tax[3,i]*PED[3,4]) + (eff_tax[4,i]*PED[4,4]))))
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
  sigma_func4 <<- approxfun(sigma_data[,c(1,5)], rule = 2)
  
  #Run the model 
  out <- data.frame(ode(y = init, func = func, times = times, parms = parms)) 
  
  n_data <- ncol(out)-1
  
  timing <- t(sapply(1:n_data, function(x)  out[max(which(!is.na(out[,x+1]))),]))
  
  if(timing[1,1] != tail(times,1)) {
    for(i in 1:n_data){
      out[seq(timing[[1]]+2, tail(times,1)+1),i+1] <- timing[i,i+1]
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
                         "R1" = data$R1 + data$R12 + data$R13 + data$R14 + data$R123 + data$R124 + data$R134 + data$R1234,
                         "R2" = data$R2 + data$R12 + data$R24 + data$R23 + data$R123 + data$R124 + data$R234 + data$R1234,
                         "R3" = data$R3 + data$R13 + data$R23 + data$R34 + data$R123 + data$R134 + data$R234 + data$R1234,
                         "R4" = data$R4 + data$R14 + data$R24 + data$R34 + data$R124 + data$R134 + data$R234 + data$R1234)
  return(agg_data)
}

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, time_between, parms, init, func, agg_func, ode_wrapper, approx_sigma){
  
  parms["time_between"] <- time_between
  
  #First Run
  run_1rd <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms, approx_sigma)[[1]])
  values_1rd <- tail(run_1rd, 1)[4:7]
  
  low_char_1rd <- names(values_1rd)[which.min(values_1rd)]
  high_char_1rd <- names(values_1rd)[which.max(values_1rd)]
  med_char_1rd <- names(values_1rd)[setdiff(1:4, c(which.min(values_1rd), which.max(values_1rd)))]
  med_char_val <- mean(as.numeric(values_1rd[med_char_1rd]))
  
  parms[["eff_tax"]][as.numeric(substr(low_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[low_char_1rd]/med_char_val)))
  parms[["eff_tax"]][as.numeric(substr(high_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[high_char_1rd]/med_char_val)))
  
  parms[["eff_tax"]][as.numeric(substr(med_char_1rd[1], 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[med_char_1rd][1]/med_char_val)))
  parms[["eff_tax"]][as.numeric(substr(med_char_1rd[2], 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[med_char_1rd][2]/med_char_val)))
  parms[["eff_tax"]][parms[["eff_tax"]] < 0.00001] <- 0
  
  #First Round of Diff Taxation
  
  parms[["int_round"]] <- 1
  
  if(int_gen > 1) {
    #All Rounds Above 1
    for(i in 2:int_gen) {
      parms[["int_round"]] <- i-1
      run <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*(i-1))), parms = parms, approx_sigma)[[1]])

      values <- tail(run, 1)[4:7]
      
      if(values[1] == 0 & values[2] == 0 & values[3] == 0) {
        parms[["eff_tax"]][,i] <- 0
      } else {
        
        low_char <- names(values)[which.min(values)]
        high_char <- names(values)[which.max(values)]
        med_char <- names(values)[setdiff(1:4, c(which.min(values), which.max(values)))]
        med_char_val <- mean(as.numeric(values[med_char]))
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/med_char_val)))

        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/med_char_val)))
        parms[["eff_tax"]][as.numeric(substr(med_char[1], 2, 2)), c((i):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[med_char][1],1)/med_char_val)))
        parms[["eff_tax"]][as.numeric(substr(med_char[2], 2, 2)), c((i):6)] <- as.numeric((parms[["base_tax"]]*(tail(run[med_char][2],1)/med_char_val)))
        
        parms[["eff_tax"]][parms[["eff_tax"]] < 0.00001] <- 0
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
  values_1rd <- tail(run_1rd, 1)[4:7]
  
  res_order_vec <- c(names(values_1rd)[which.max(values_1rd)],
                     names(values_1rd)[setdiff(1:4, c(which.min(values_1rd), which.max(values_1rd)))][1],
                     names(values_1rd)[setdiff(1:4, c(which.min(values_1rd), which.max(values_1rd)))][2],
                     names(values_1rd)[which.min(values_1rd)])[res_order]
  
  parms[["eff_tax"]][as.numeric(substr(res_order_vec, 2, 2)), c(1:6)] <- as.numeric(parms[["base_tax"]])
  parms[["int_round"]] <- 1
  
  #Real Model Run 
  run_real <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms, approx_sigma)
  return(run_real)
}

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, 
          R1 = 0, R2 = 0, R3 = 0, R4 = 0,
          R12 = 0, R13 = 0, R14 = 0, R23 = 0, R24 = 0, R34 = 0,
          R123 = 0, R124 = 0, R134 = 0, R234 = 0, 
          R1234 = 0)

parms = list(lambda = 1/365*(2), int_round = 1, 
             beta = 5, sigma1 = 0.2, sigma2 = 0.2, sigma3 = 0.2, sigma4 = 0.2,
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, r_rrrr = 1/7, 
             r_t = 1/6, eta_wr = 0.3, eta_rw = 0.04, 
             eta_rr = 0.01, eta_rrr = 0.01, eta_rrrr = 0.01,  
             c1 = 0.945, c2 = 0.92, c3 = 0.865, c4 = 0.8,
             c12 = 0.8, c13 = 0.76,  c14 = 0.75, 
             c23 = 0.75, c24 = 0.71, c34 = 0.7,
             c123 = 0.7, c124 = 0.65, c134 = 0.625, c234 = 0.6,
             c1234 = 0.6,
             PED = matrix(c(-1, 0.25, 0.25, 0.25, 
                            0.25, -1, 0.25, 0.25,
                            0.25, 0.25, -1, 0.25,
                            0.25, 0.25, 0.25, -1), #Be aware of this matrix
                          nrow = 4, ncol = 4, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0), 
                              nrow = 4, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# Baseline Model ----------------------------------------------------------

parms1 <- parms
testrun_flat <- ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms1, approx_sigma)[[1]]
test_run_agg <- agg_func(testrun_flat) 

test_run_agg$AverageRes <- rowMeans(test_run_agg[,4:7])
test_run_agg$TotInf <- rowSums(testrun_flat[,3:18])

parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5; parms1[["int_round"]] <- 1
testrun_flat <- ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms1, approx_sigma)[[1]]
ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms1, approx_sigma)[[2]]
test_plot_flat <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[4:7])
ggplot(test_plot_flat, aes(time, value, color = variable)) + geom_line() + theme_bw()

# Intervention Scenarios --------------------------------------------------

#Flat Taxs
parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5; parms1[["int_round"]] <- 1
testrun_flat <- list(ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms1, approx_sigma)[[1]])

#Single Tax 
single_list <- list()

for(i in 1:4) {
  parms1 <- parms
  single_list[[i]] <- single_tax(i, 0.5, parms1, init, amr, agg_func, ode_wrapper, approx_sigma)[[1]]
  single_tax(i, 0.5, parms1, init, amr, agg_func, ode_wrapper, approx_sigma)[[2]]
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  diff_tax_list[[i]] <- multi_int_fun(i, 365*3, parms, init, amr, agg_func, ode_wrapper, approx_sigma)[[1]]
}

# Plotting the Scenarios --------------------------------------------------

#Create a combined list of all the scenarios
list_scen <- unlist(list(testrun_flat, single_list, diff_tax_list), recursive = FALSE)

melt_data <- list()

#Melt each one
for(i in 1:length(list_scen)) {
  data_agg <- agg_func(list_scen[[i]]) 
  colnames(data_agg)[4:7] <- c("High Res (HR)", "Medium Res (HR)", "Medium Res (MR2)",  "Low Res (LR)") 
  melt_data[[i]] <- data.frame(melt(data_agg, id.vars = "time", measure.vars = colnames(data_agg)[4:7]),
                               "scen" = c("flat", "singleHR", "singleMR1","singleMR2","singleLR",
                                          "diff1", "diff2", "diff3", "diff4", "diff5", "diff6")[i])
}

p_data <- list()

#Plotting Loop
for(i in 1:length(melt_data)) {
  data <- melt_data[[i]]
  p_data[[i]] <- ggplot(data, aes(time, value, color = variable)) + geom_line() + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "Time", y = "Prevalence", title = c("Flat Tax",
                                                 "Single Tax (HR)","Single Tax (MR1)", "Single Tax (MR2)", "Single Tax (LR)",
                                                 "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)",
                                                 "Diff Tax (4 Rd)", "Diff Tax (5 Rd)", "Diff Tax (6 Rd)")[i], color = "")
}

ggarrange(p_data[[1]], "", "", "",
          p_data[[2]], p_data[[3]], p_data[[4]], p_data[[5]],
          p_data[[6]], p_data[[7]], p_data[[8]], p_data[[9]],
          p_data[[10]], p_data[[11]], "" , "",
          labels = c("A", "", "", "",
                     "B", "", "", "",
                     "C", "", "", "",
                      "", "", "", ""), hjust = -.1,  nrow = 4, ncol = 4, common.legend = T, legend = "bottom")

library("deSolve"); library("parallel")
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
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i):6)] <- (((1 + as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/med_char_val)))) - 
                                                                                (1 + parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)])) / 
                                                                               (1 + parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)]))
        
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i):6)] <- (((1 + as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/med_char_val)))) - 
                                                                                 (1 + parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)])) / 
                                                                                (1 + parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)]))
        
        parms[["eff_tax"]][as.numeric(substr(med_char[1], 2, 2)), c((i):6)] <- (((1 + as.numeric((parms[["base_tax"]]*(tail(run[med_char][1],1)/med_char_val)))) - 
                                                                                   (1 + parms[["eff_tax"]][as.numeric(substr(med_char[1], 2, 2)), (i-1)])) / 
                                                                                  (1 + parms[["eff_tax"]][as.numeric(substr(med_char[1], 2, 2)), (i-1)]))
        
        parms[["eff_tax"]][as.numeric(substr(med_char[2], 2, 2)), c((i):6)] <- (((1 + as.numeric((parms[["base_tax"]]*(tail(run[med_char][2],1)/med_char_val)))) - 
                                                                                   (1 + parms[["eff_tax"]][as.numeric(substr(med_char[2], 2, 2)), (i-1)])) / 
                                                                                  (1 + parms[["eff_tax"]][as.numeric(substr(med_char[2], 2, 2)), (i-1)]))
      }
      parms[["int_round"]] <- i
    }
  }
  out_run <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms, approx_sigma)
  return(out_run)
}

# Aggregated Function -----------------------------------------------------

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

# Extract Usage -----------------------------------------------------------

usage_fun <- function(parms){
  
  #In case there is no intervention
  if(parms[["int_round"]] == 0) {
    usage = data.frame("time" = seq(0,7000),
                       "PopUsage1" = rep(parms[["sigma_mat"]][1,2], 7001),
                       
                       "PopUsage2" = rep(parms[["sigma_mat"]][2,2], 7001),
                       
                       "PopUsage3" = rep(parms[["sigma_mat"]][3,2], 7001),
                       
                       "PopUsage4" = rep(parms[["sigma_mat"]][4,2], 7001))
  }
  
  #Intervention
  if(parms[["int_round"]] > 0) {
    usage = data.frame("time" = seq(0,7000),
                       "PopUsage1" = c(rep(parms[["sigma_mat"]][1,2] , 365*3), rep(parms[["sigma_mat"]][1,3] , 365*3), rep(parms[["sigma_mat"]][1,4] , 365*3),
                                       rep(parms[["sigma_mat"]][1,5], 365*3), rep(parms[["sigma_mat"]][1,6] , 365*3), rep(parms[["sigma_mat"]][1,7] , (365*3)+(7001-(365*3)*6))),
                       
                       "PopUsage2" = c(rep(parms[["sigma_mat"]][2,2] , 365*3), rep(parms[["sigma_mat"]][2,3] , 365*3), rep(parms[["sigma_mat"]][2,4] , 365*3),
                                       rep(parms[["sigma_mat"]][2,5], 365*3), rep(parms[["sigma_mat"]][2,6] , 365*3), rep(parms[["sigma_mat"]][2,7] , (365*3)+(7001-(365*3)*6))),
                       
                       "PopUsage3" = c(rep(parms[["sigma_mat"]][3,2] , 365*3), rep(parms[["sigma_mat"]][3,3] , 365*3), rep(parms[["sigma_mat"]][3,4] , 365*3),
                                       rep(parms[["sigma_mat"]][3,5], 365*3), rep(parms[["sigma_mat"]][3,6] , 365*3), rep(parms[["sigma_mat"]][3,7] , (365*3)+(7001-(365*3)*6))),
                       
                       "PopUsage4" = c(rep(parms[["sigma_mat"]][4,2] , 365*3), rep(parms[["sigma_mat"]][4,3] , 365*3), rep(parms[["sigma_mat"]][4,4] , 365*3),
                                       rep(parms[["sigma_mat"]][4,5], 365*3), rep(parms[["sigma_mat"]][4,6] , 365*3), rep(parms[["sigma_mat"]][4,7] , (365*3)+(7001-(365*3)*6))))
  }
  
  usage$totusage = rowSums(usage[2:5])
  usage$reduc_use = parms[["sigma1"]] + parms[["sigma2"]] + parms[["sigma3"]] + parms[["sigma4"]]  - usage$totusage
  
  return(usage)
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
             c1 = 0.945, c2 = 0.935, c3 = 0.91, c4 = 0.8,
             c12 = 0.8, c13 = 0.775,  c14 = 0.75, 
             c23 = 0.75, c24 = 0.725, c34 = 0.7,
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

# The Function ------------------------------------------------------------

low_parm <- c(1/3650*(2), #lambda
              0, #beta
              0, #sigma1
              0, #sigma2
              0, #sigma3
              0, #sigma4
              1/50, #r_wt
              1/50, #r_r
              1/50, #r_rr
              1/50, #r_rrr
              1/50, #r_rrrr
              1/50, #r_t
              0.03, #eta_wr
              0.004, #eta_rw
              0.001, #eta_rr
              0.001, #eta_rrr
              0.001, #eta_rrrr
              0.5, #c1
              0.5, #c2
              0.5, #c3
              0.5, #c4
              0.5, #c12
              0.5, #c13
              0.5, #c14
              0.5, #c23
              0.5, #c24
              0.5, #c34
              0.5, #c123
              0.5, #c124
              0.5, #c134
              0.5, #c234
              0.5, #c1234
              0, #rho
              0) #baseline tax

high_parm <- c(1/36.5*(2), #lambda
               10, #beta
               1, #sigma1
               1, #sigma2
               1, #sigma3
               1, #sigma4
               1/2, #r_wt
               1/2, #r_r
               1/2, #r_rr
               1/2, #r_rrr
               1/2, #r_rrrr
               1/2, #r_t
               3, #eta_wr
               0.4, #eta_rw
               0.1, #eta_rr
               0.1, #eta_rrr
               0.1, #eta_rrrr
               1, #c1
               1, #c2
               1, #c3
               1, #c4
               1, #c12
               1, #c13
               1, #c14
               1, #c23
               1, #c24
               1, #c34
               1, #c123
               1, #c124
               1, #c134
               1, #c234
               1, #c1234
               1, #rho
               1) #baseline tax

#Creating the Parm Dataframe

parm_data <- data.frame(t(replicate(10000, runif(34, low_parm, high_parm))))

colnames(parm_data) <- c("lambda", "beta", "sigma1", "sigma2", "sigma3",  "sigma4", 
                         "r_wt", "r_r", "r_rr", "r_rrr","r_rrrr","r_t",
                         "eta_wr", "eta_rw", "eta_rr", "eta_rrr", "eta_rrrr",
                         "c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                         "c123", "c124","c134", "c234", "c1234",  
                         "rho", "base_tax")

for(i in 1:nrow(parm_data)) {
  if(sum(parm_data[c("sigma1", "sigma2", "sigma3", "sigma4")][i,]) > 1) {
    parm_data[c("sigma1", "sigma2", "sigma3", "sigma4")][i,] <- parm_data[c("sigma1", "sigma2", "sigma3", "sigma4")][i,]/
      (sum(parm_data[c("sigma1", "sigma2", "sigma3", "sigma4")][i,]) + runif(1, 0, 1))
  }
}

parm_data[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr", "eta_rrrr")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr", "eta_rrrr")][x,]), decreasing = T)))

parm_data[c("r_wt", "r_r", "r_rr", "r_rrr", "r_rrrr", "r_t")] <- t(sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("r_wt", "r_r", "r_rr", "r_rrr", "r_rrrr", "r_t")][x,]), decreasing = F)))


parm_data[c("c1", "c2", "c3", "c4")] <- t(sapply(1:nrow(parm_data), function(x) 
  sample(sort(as.numeric(parm_data[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                                     "c123", "c124","c134", "c234", "c1234")][x,]), decreasing = T)[1:4], 
         size = 4, replace = FALSE)))

parm_data[c("c12", "c13", "c14", "c23","c24","c34")] <-  t(sapply(1:nrow(parm_data), function(x) 
  sample(sort(as.numeric(parm_data[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                                     "c123", "c124","c134", "c234", "c1234")][x,]), decreasing = T)[5:10], 
         size = 6, replace = FALSE)))

parm_data[c("c123", "c124","c134", "c234")] <-  t(sapply(1:nrow(parm_data), function(x) 
  sample(sort(as.numeric(parm_data[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                                     "c123", "c124","c134", "c234", "c1234")][x,]), decreasing = T)[11:14], 
         size = 4, replace = FALSE)))

parm_data["c1234"] <- sapply(1:nrow(parm_data), function(x) 
  sort(as.numeric(parm_data[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                              "c123", "c124","c134", "c234", "c1234")][x,]), decreasing = T)[15])

parm_data_comb <- data.frame(parm_data, t_n = 3000, int_round = 0,
                             time_between = Inf)

# Creating the Parallel Montonicity Function ------------------------------

mono_func <- function(n, parms_frame, init, amr_ode, usage_fun, multi_int_fun, low_parm, high_parm, agg_func, thresh, approx_sigma) {
  
  parms_base = as.list(parms_frame[n,])
  parms_base = append(parms_base, parms["PED"])
  parms_base = append(parms_base, parms["eff_tax"])
  
  #Run Baseline
  run_base <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, approx_sigma)[[1]]
  run_base_agg <- agg_func(run_base)
  values <- tail(run_base_agg, 1)
  
  if(values[4] == 0 & values[5] == 0 & values[6] == 0 & values[7] == 0) {
    while(values[4] == 0 & values[5] == 0 & values[6] == 0 & values[7] == 0) {
      parms_base[c(1:34)] <- as.list(runif(34, low_parm, high_parm))
      
      if(sum(unlist(parms_base[c("sigma1", "sigma2", "sigma3", "sigma4")])) > 1) {
        parms_base[c("sigma1", "sigma2", "sigma3", "sigma4")] <- as.list(unlist(parms_base[c("sigma1", "sigma2", "sigma3", "sigma4")])/
                                                                 (sum(unlist(parms_base[c("sigma1", "sigma2", "sigma3",  "sigma4")])) + runif(1, 0, 1)))
      }
      
      parms_base[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr", "eta_rrrr")] <- 
        as.list(sort(as.numeric(parms_base[c("eta_wr", "eta_rw", "eta_rr", "eta_rrr", "eta_rrrr")]), decreasing = T))
      
      parms_base[c("r_wt", "r_r", "r_rr", "r_rrr","r_rrrr", "r_t")] <- 
        as.list(sort(as.numeric(parms_base[c("r_wt", "r_r", "r_rr", "r_rrr","r_rrrr", "r_t")]), decreasing = F))
      
      parms_base[c("c1", "c2", "c3", "c4")] <- 
        as.list(sample(sort(as.numeric(parms_base[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                                                    "c123", "c124","c134", "c234", "c1234")]), decreasing = T)[1:4]), size = 4, replace = F)
      
      parms_base[c("c12", "c13", "c14", "c23","c24","c34")] <- 
        as.list(sample(sort(as.numeric(parms_base[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                                                    "c123", "c124","c134", "c234", "c1234")]), decreasing = T)[5:10]), size = 6, replace = F)
      
      parms_base[c("c123", "c124","c134", "c234")] <- 
        as.list(sample(sort(as.numeric(parms_base[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                                                    "c123", "c124","c134", "c234", "c1234")]), decreasing = T)[11:14]), size = 4, replace = F)
      
      parms_base["c1234"] <- 
        as.list(sort(as.numeric(parms_base[c("c1", "c2", "c3","c4", "c12", "c13", "c14", "c23","c24","c34", 
                                             "c123", "c124","c134", "c234", "c1234")]), decreasing = T)[15])
      
      run_base <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10000), parms = parms_base, approx_sigma)[[1]]
      run_base_agg <- agg_func(run_base)
      values <- tail(run_base_agg, 1)
    }
  }

  run <- run_base[run_base[,1] > parms_base[["t_n"]],]
  run_base_agg <- run_base_agg[run_base_agg[,1] > parms_base[["t_n"]],]
  
  #Identifying the order of the resistances
  res_order_vec <- c(names(values[4:7])[which.max(values[4:7])],
                     names(values[4:7])[setdiff(1:4, c(which.min(values[4:7]), which.max(values[4:7])))][1],
                     names(values[4:7])[setdiff(1:4, c(which.min(values[4:7]), which.max(values[4:7])))][2],
                     names(values[4:7])[which.min(values[4:7])])
  
  base_tot_inf <- signif(sum(run[3:18]), 5)
  base_int_res <- signif(sum(rowMeans(run_base_agg[4:7])), 5)
  
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
    if(i >= 2 & i <= 5) {
      parms[["eff_tax"]][as.numeric(substr(res_order_vec[i-1], 2, 2)), c(1:6)] <- parms[["base_tax"]]
      parms[["int_round"]] <- 1
      out_run <- ode_wrapper(y = init, func = amr_ode, times = seq(0, 10000), parms = parms, approx_sigma)
      out <- out_run[[1]]
      parms <- out_run[[2]]
    }
    if(i >= 6 & i <= 11) {
      diff <- multi_int_fun(i-5, 365*3, parms, init, amr_ode, agg_func, ode_wrapper, approx_sigma)
      out <- diff[[1]]
      parms <- diff[[2]]
    }
    
    data_temp <- out[out[,1] > parms[["t_n"]],]
    data_temp_agg <- agg_func(data_temp)
    
    out_vec <- signif(c(sum(data_temp[3:18]),
                        sum(rowMeans(data_temp_agg[4:7]))),5)
    reduc_usage_vec <- sum(usage_fun(parms)[,7])
    
    #Aggregation
    
    out$aggR1 <- out$R1 + out$R12 + out$R13 + out$R14 + out$R123 + out$R124 + out$R134 + out$R1234
    out$aggR2 <- out$R2 + out$R12 + out$R24 + out$R23 + out$R123 + out$R124 + out$R234 + out$R1234
    out$aggR3 <- out$R3 + out$R13 + out$R23 + out$R34 + out$R123 + out$R134 + out$R234 + out$R1234
    out$aggR4 <- out$R4 + out$R14 + out$R24 + out$R34 + out$R124 + out$R134 + out$R234 + out$R1234
    
    #Determine the X% Thresholds that you want to be under
    thresholds <- unlist(out[parms[["t_n"]]-1, 19:22]*thresh)
    under_thresh <- sweep(out[out[,1] > parms[["t_n"]],][,19:22], 2, thresholds)
    
    #Calculate the number of days you are under said threshold
    under_50 <- c(nrow(under_thresh[under_thresh$aggR1 < 0,]), 
                  nrow(under_thresh[under_thresh$aggR2 < 0,]), 
                  nrow(under_thresh[under_thresh$aggR3 < 0,]),
                  nrow(under_thresh[under_thresh$aggR4 < 0,]))

    #Find the Sum and make each value proportionate to one another 

    prop_vec <- sum(under_50) / (10000 - parms[["t_n"]])
    
    prop_vec_shan <- under_50 / sum(under_50)
    prop_vec_shan <- prop_vec_shan[prop_vec_shan != 0]
    
    #Store Computation Vectors 
    store_vec_inf[i] <- (out_vec[1] - base_tot_inf)/reduc_usage_vec
    store_vec_res[i] <- (base_int_res - out_vec[2])/reduc_usage_vec
    store_vec_shan[i] <- -sum(sapply(1:length(prop_vec_shan), function(x) prop_vec_shan[x]*log(prop_vec_shan[x])))
    store_vec_avganti[i] <- prop_vec
  }
  
  output <- c(store_vec_inf, store_vec_res, store_vec_shan,store_vec_avganti, parms_base[c(1:38)])
  names(output) <- c("flat_inf", "singleHR_inf", "singleMR1_inf", "singleMR2_inf", "singleLR_inf", "diff1_inf", "diff2_inf", "diff3_inf", "diff4_inf", "diff5_inf", "diff6_inf", 
                     "flat_res", "singleHR_res", "singleMR1_res", "singleMR2_res", "singleLR_res", "diff1_res", "diff2_res", "diff3_res", "diff4_res", "diff5_res", "diff6_res",
                     "flat_shan", "singleHR_shan", "singleMR1_shan", "singleMR2_shan", "singleLR_shan", "diff1_shan", "diff2_shan", "diff3_shan", "diff4_shan", "diff5_shan", "diff6_shan",
                     "flat_avganti", "singleHR_avganti", "singleMR_avganti", "singleMR2_avganti", "singleLR_avganti", "diff1_avganti", "diff2_avganti", "diff3_avganti", "diff4_avganti", "diff5_avganti", "diff6_avganti",
                     names(parms_base[c(1:38)]))
  
  return(output)
}

# Run the Model ----------------------------------------------------------

start_time <- Sys.time()

test <- mclapply(1:1000, 
                 FUN = mono_func, 
                 parms_frame = parm_data_comb, 
                 init = c(X = 0.99, Wt = 1-0.99, 
                          R1 = 0, R2 = 0, R3 = 0, R4 = 0,
                          R12 = 0, R13 = 0, R14 = 0, R23 = 0, R24 = 0, R34 = 0,
                          R123 = 0, R124 = 0, R134 = 0, R234 = 0, 
                          R1234 = 0), 
                 amr_ode = amr, 
                 usage_fun = usage_fun,
                 multi_int_fun = multi_int_fun,
                 low_parm = low_parm,
                 high_parm = high_parm,
                 agg_func = agg_func,
                 thresh = 0.5,
                 mc.cores = 10, 
                 approx_sigma = approx_sigma) 

#Combine the Output into a "normal" looking dataframe
comb_data <- data.frame(do.call(rbind, test))
comb_data_new <- data.frame(matrix(NA, nrow = nrow(comb_data), ncol = 44))

for(i in 1:nrow(comb_data)) {
  comb_data_new[i,] <- unlist(comb_data[i,1:44])
}

colnames(comb_data_new) <- colnames(comb_data)[1:44]

#Update the Parameter Set 
parm_data_comb_new <- parm_data_comb
parm_data_comb_new[1:nrow(comb_data),1:34] <- comb_data[,45:78]

parm_list <- list()

for(i in 1:nrow(parm_data_comb_new)) {
  p_list <- as.list(unlist(parm_data_comb_new[i,]))
  p_list <- append(p_list, parms["eff_tax"])
  p_list <- append(p_list, parms["PED"])
  parm_list[[i]] <- p_list
}

#Save the output

saveRDS(parm_list, "/cluster/home/amorgan/Sens_Anal_Output/MDR_run_parms_four.RDS")
saveRDS(comb_data_new, "/cluster/home/amorgan/Sens_Anal_Output/MDR_run_four.RDS")

end_time <- Sys.time()
print(end_time - start_time)

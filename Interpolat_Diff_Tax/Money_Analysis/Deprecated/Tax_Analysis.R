library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr")
rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Model_Fit/Model_Output/New")

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
        
        parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), c((i):6)] <- (((1+as.numeric((parms[["base_tax"]]*(tail(run[low_char],1)/values_1rd[med_char_1rd]))))-
                                                                                (1+parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)]))/
                                                                               (1+parms[["eff_tax"]][as.numeric(substr(low_char, 2, 2)), (i-1)]))
        parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), c((i):6)] <- (((1+as.numeric((parms[["base_tax"]]*(tail(run[high_char],1)/values_1rd[med_char_1rd]))))-
                                                                                 (1+parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)]))/
                                                                                (1+parms[["eff_tax"]][as.numeric(substr(high_char, 2, 2)), (i-1)]))
        parms[["eff_tax"]][as.numeric(substr(med_char, 2, 2)), c((i):6)] <- (((1+as.numeric((parms[["base_tax"]]*(tail(run[med_char],1)/values_1rd[med_char_1rd]))))-
                                                                                (1+parms[["eff_tax"]][as.numeric(substr(med_char, 2, 2)), (i-1)]))/
                                                                               (1+parms[["eff_tax"]][as.numeric(substr(med_char, 2, 2)), (i-1)]))
      }
      parms[["int_round"]] <- i
    }
  }
  out_run <- ode_wrapper(y = init, func = func, times = seq(0, 10300), parms = parms, approx_sigma)
  return(out_run)
}

# Baseline Parms ----------------------------------------------------------

post_dist_names <- grep("ABC_v3_",
                        list.files("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Model_Fit/Model_Output/New"), value = TRUE)

post_dist <- lapply(post_dist_names, read.csv)

post_dist <- mapply(cbind, post_dist, "gen" = sapply(1:length(post_dist), function(x) paste0("gen_", x)), 
                    SIMPLIFY=F)
post_dist <- do.call("rbind", post_dist)

maps_est <- colMeans(post_dist[post_dist$gen == tail(unique(post_dist$gen),1),][,1:11])

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = list(lambda = 1/365*(2), int_round = 1, 
             beta = maps_est["beta"], 
             sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
             r_t = 1/7, 
             eta_wr = maps_est["eta_wr"], 
             eta_rw = maps_est["eta_rw"], 
             eta_rr = maps_est["eta_rr_rrr"], eta_rrr = maps_est["eta_rr_rrr"],  
             c1 = maps_est["c1"], c2 = maps_est["c2"], 
             c3 = maps_est["c3"],
             c12 = maps_est["c12"], c13 = maps_est["c13"], 
             c23 = maps_est["c23"],
             c123 = maps_est["c123"],
             PED = matrix(c(-1.5, 0.75, 0.5, 
                            0.25, -1.25, 0.75,
                            0, 0.25, -1), #Be aware of this matrix
                          nrow = 3, ncol = 3, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 3, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# Intervention Scenarios --------------------------------------------------

#Flat Tax

parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5; parms1[["int_round"]] <- 1
testrun_flat <- list(ode_wrapper(y = init, func = amr, times = seq(0, 10300), parms = parms1, approx_sigma)[[2]])

#Single Tax 
single_list <- list()

for(i in 1:3) {
  parms1 <- parms
  single_list[[i]] <- single_tax(i, 0.5, parms1, init, amr, agg_func, ode_wrapper, approx_sigma)[[2]]
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  diff_tax_list[[i]] <- multi_int_fun(i, 365*3, parms, init, amr, agg_func, ode_wrapper, approx_sigma)[[2]]
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

#Testing the Function
efftax_mat <- diff_tax_list[[3]]$eff_tax
test_mat <- conv_diff_efftax(efftax_mat,3)

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

tax_data <- year_tax(test_mat)

# Drug Class Dataframe ----------------------------------------------------

data_anti <- data.frame("Antibiotic" = c("Tetracyclines", "Amphenicols", "Penicillin", "Cephalosporins", "Sulphonamides", "Macrolides",
                                         "Aminoglycosides", "Quinolones", "Lincosamides"),
                        "Revenue" = c(4391696719, 132574205.9, 331097866.7, 512035885, 368451522.6, 10144501366, 4609169528,129688300.4, 136177433.3),
                        "Resistance" = c(31.57, 6.9, 9.13, 6.10, 19.80, 0.73, 9.10, 3.53, 0),
                        "Group" = c(1, 2, 2, 2, 1, 3, 2, 3, 3))

tax_vector <- c("G1" = sum(data_anti[data_anti$Group == 1,]$Revenue), 
                "G2" = sum(data_anti[data_anti$Group == 2,]$Revenue), 
                "G3" = sum(data_anti[data_anti$Group == 3,]$Revenue))

# Identifying Tax ---------------------------------------------------------

list_scen <- unlist(list(testrun_flat, single_list, diff_tax_list), recursive = FALSE)
tax_list <- list()

for(i in 1:10){
  eff_tax <- list_scen[[i]]$eff_tax
  
  if(i > 5) {
    eff_tax <- conv_diff_efftax(list_scen[[i]]$eff_tax, i-4)
  }
  
  data_year <- year_tax(eff_tax)
  
  data_year_new <- data_year
  data_year_new[,2:4] <- sweep(data_year_new[,2:4], 2, tax_vector, "*")
  data_year_new$total <- rowSums(data_year_new[,2:4])
  
  tax_list[[i]] <- data_year_new
}


library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
library("rbenchmark")

rm(list=ls())

# ODEs --------------------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    #Set up the matrix for the Sigmas
    
    sigma_use1 <- sigma_mat[1,1]; sigma_use2 <- sigma_mat[2,1]; sigma_use3 <- sigma_mat[3,1]  
    
    #For One Round
    if(t > t_n & int_round > 0) { 
      sigma_use1 <- sigma_mat[1,2]; sigma_use2 <- sigma_mat[2,2]; sigma_use3 <- sigma_mat[3,2]  
    }
    
    if(t > (t_n + time_between) & int_round > 1) { 
      sigma_use1 <- sigma_mat[1,3]; sigma_use2 <- sigma_mat[2,3]; sigma_use3 <- sigma_mat[3,3]  
    }
    
    if(t > (t_n + time_between*2) & int_round > 2) { 
      sigma_use1 <- sigma_mat[1,4]; sigma_use2 <- sigma_mat[2,4]; sigma_use3 <- sigma_mat[3,4]  
    }
    
    if(t > (t_n + time_between*3) & int_round > 3) {
      sigma_use1 <- sigma_mat[1,5]; sigma_use2 <- sigma_mat[2,5]; sigma_use3 <- sigma_mat[3,5]  
    }
    
    if(t > (t_n + time_between*4) & int_round > 4) { 
      sigma_use1 <- sigma_mat[1,6]; sigma_use2 <- sigma_mat[2,6]; sigma_use3 <- sigma_mat[3,6]  
    }
    
    if(t > (t_n + time_between*5) & int_round > 5) {
      sigma_use1 <- sigma_mat[1,7]; sigma_use2 <- sigma_mat[2,7]; sigma_use3 <- sigma_mat[3,7]  
    }
    
    #ODES Below
    
    dX = lambda - lambda*X - beta*X*(Wt + R1*c1 + R2*c2 + R3*c3 + R12*c12 + R13*c13 + R23*c23 + R123*c123) +
      r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + r_r*R1*sigma_use1 + r_r*R2*sigma_use2 + r_r*R3*sigma_use3 +
      r_rr*R12*(sigma_use1 + sigma_use2) + r_rr*R13*(sigma_use1 + sigma_use3) + r_rr*R23*(sigma_use2 + sigma_use3) + 
      r_rrr*R123*(sigma_use1 + sigma_use2 + sigma_use3) + r_t*(1-rho)*(Wt*(sigma_use1 + sigma_use2 + sigma_use3) + R1*(sigma_use2 + sigma_use3) + 
                                                                         R2*(sigma_use1 + sigma_use3) + R3*(sigma_use1 + sigma_use2) + 
                                                                         R12*sigma_use3 + R13*sigma_use2 + R23*sigma_use1)
    
    dWt = - lambda*Wt + beta*X*Wt - r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) - r_t*Wt*(1-rho)*(sigma_use1 + sigma_use2 + sigma_use3) +
      eta_rw*(R1 + R2 + R3 + R12 + R13 + R23 + R123)*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) - 
      eta_wr*Wt*rho*(sigma_use1 + sigma_use2 + sigma_use3)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma_use2 + sigma_use3)*R1 - r_r*sigma_use1*R1 - eta_rr*R1*rho*sigma_use2 -
      eta_rr*R1*rho*sigma_use3 - eta_rw*R1*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma_use1 + sigma_use3)*R2 - r_r*sigma_use2*R2 - eta_rr*R2*rho*sigma_use1 -
      eta_rr*R2*rho*sigma_use3 - eta_rw*R2*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use2
    
    dR3 = - lambda*R3 + beta*X*R3*c3 - r_t*(1-rho)*(sigma_use1 + sigma_use2)*R3 - r_r*sigma_use3*R3 - eta_rr*R3*rho*sigma_use1 -
      eta_rr*R3*rho*sigma_use2 - eta_rw*R3*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use3
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_t*(1-rho)*sigma_use3*R12 - r_rr*(sigma_use1 + sigma_use2)*R12 - eta_rrr*R12*rho*sigma_use3 -
      eta_rw*R12*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R1*rho*sigma_use2 + eta_rr*R2*rho*sigma_use1 
    
    dR13 = - lambda*R13 + beta*X*R13*c13 - r_t*(1-rho)*sigma_use2*R13 - r_rr*(sigma_use1 + sigma_use3)*R13 - eta_rrr*R13*rho*sigma_use2 -
      eta_rw*R13*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R1*rho*sigma_use3 + eta_rr*R3*rho*sigma_use1 
    
    dR23 = - lambda*R23 + beta*X*R23*c23 - r_t*(1-rho)*sigma_use1*R23 - r_rr*(sigma_use2 + sigma_use3)*R23 - eta_rrr*R23*rho*sigma_use1 -
      eta_rw*R23*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R2*rho*sigma_use3 + eta_rr*R3*rho*sigma_use2 
    
    dR123 = - lambda*R123 + beta*X*R123*c123 - r_rrr*(sigma_use1 + sigma_use2 + sigma_use3)*R123 - eta_rw*R123*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + 
      eta_rrr*rho*(sigma_use3*R12 + sigma_use2*R13 + sigma_use1*R23)
    
    return(list(c(dX,dWt,
                  dR1,dR2,dR3,
                  dR12,dR13,dR23,
                  dR123)))
  })
}

# ODE Function Wrapper ----------------------------------------------------

ode_wrapper <- function(times, y, parms, func) {
  
  sigma_mat = matrix(c(rep(parms[["sigma1"]], 7),
                       rep(parms[["sigma2"]], 7),
                       rep(parms[["sigma3"]], 7)), 
                     nrow = 3, ncol = 7, byrow = T)
  
  eff_tax <- parms[["eff_tax"]]; PED <- parms[["PED"]]
  
  if(parms[["int_round"]] > 0 ) {
    for(i in 1:parms[["int_round"]]) {
      stor_sigma <- sigma_mat[,i]
      if(sigma_mat[1,i] == 0 | sigma_mat[2,i] == 0 | sigma_mat[3,i] == 0) {
        stor_sigma[stor_sigma == 0] = 0.01
      }
      
      sigma_mat[,(i+1):7] = c(stor_sigma[1]*(1 + ((eff_tax[1,i]*PED[1,1]) + (eff_tax[2,i]*PED[2,1]) + (eff_tax[3,i]*PED[3,1]))),
                              stor_sigma[2]*(1 + ((eff_tax[1,i]*PED[1,2]) + (eff_tax[2,i]*PED[2,2]) + (eff_tax[3,i]*PED[3,2]))),
                              stor_sigma[3]*(1 + ((eff_tax[1,i]*PED[1,3]) + (eff_tax[2,i]*PED[2,3]) + (eff_tax[3,i]*PED[3,3]))))
      sigma_mat[sigma_mat < 0] = 0
      
      if(colSums(sigma_mat)[i+1] > 1) {
        sigma_mat[,(i+1):7] <- sigma_mat[,i+1]/colSums(sigma_mat)[i+1]
      }
    }
  }
  
  parms[["sigma_mat"]] <- sigma_mat
  
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

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, time_between, parms, init, func, agg_func, ode_wrapper){
  
  parms["time_between"] <- time_between
  
  #First Run
  run_1rd <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms)[[1]])
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  low_char_1rd <- names(values_1rd)[which.min(values_1rd)]
  high_char_1rd <- names(values_1rd)[which.max(values_1rd)]
  med_char_1rd <- names(values_1rd)[setdiff(1:3, c(which.min(values_1rd), which.max(values_1rd)))]
  
  parms[["eff_tax"]][as.numeric(substr(low_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[low_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][as.numeric(substr(high_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[high_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][as.numeric(substr(med_char_1rd, 2, 2)), c(1:6)] <- as.numeric((parms[["base_tax"]]*(values_1rd[med_char_1rd]/values_1rd[med_char_1rd])))
  parms[["eff_tax"]][parms[["eff_tax"]] < 0] <- 0
  
  parms[["int_round"]] <- 1
  
  if(int_gen > 1) {
    
    for(i in 1:(int_gen-1)) {
      parms[["int_round"]] <- int_gen
      run <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*i)), parms = parms)[[1]])
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
  out_run <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms)
  return(out_run)
}

# Single Taxation Function ------------------------------------------------

single_tax <- function(res_order, tax, parms, init, func, agg_func, ode_wrapper) {
  
  #First Run
  parms[["base_tax"]] <- tax
  
  run_1rd <- agg_func(ode_wrapper(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms)[[1]])
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  res_order_vec <- c(names(values_1rd)[which.max(values_1rd)],
                     names(values_1rd)[setdiff(1:3, c(which.min(values_1rd), which.max(values_1rd)))],
                     names(values_1rd)[which.min(values_1rd)])[res_order]
  
  parms[["eff_tax"]][as.numeric(substr(res_order_vec, 2, 2)), c(1:6)] <- as.numeric(parms[["base_tax"]])
  parms[["int_round"]] <- 1
  
  #Real Model Run 
  run_real <- ode_wrapper(y = init, func = func, times = seq(0, 10000), parms = parms)
  return(run_real)
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


# Extract Usage -----------------------------------------------------------

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

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
          R12 = 0, R13 = 0, R23 = 0,
          R123 = 0)

parms = list(lambda = 1/365*(2), int_round = 1, 
             beta = 5, sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
             r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
             r_t = 1/7, eta_wr = 0.3, eta_rw = 0.04, 
             eta_rr = 0.01, eta_rrr = 0.01,  
             c1 = 0.945, c2 = 0.925, c3 = 0.85,
             c12 = 0.845, c13 = 0.825, c23 = 0.75,
             c123 = 0.7,
             PED = matrix(c(-1, 0.4, 0.4, 
                            0.4, -1, 0.4,
                            0.4, 0.4, -1), #Be aware of this matrix
                          nrow = 3, ncol = 3, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 3, ncol = 6, byrow = T),
             sigma_mat = matrix(c(0, 0, 0, 0, 0, 0, 0, 
                                  0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0), 
                                nrow = 3, ncol = 7, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# Intervention Scenarios --------------------------------------------------

#Flat Tax
parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5; parms1[["int_round"]] <- 1
testrun_flat <- list(ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms1))

#Single Tax 
single_list <- list()

for(i in 1:3) {
  parms1 <- parms
  single_list[[i]] <- single_tax(i, 0.5, parms1, init, amr, agg_func, ode_wrapper)
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  diff_tax_list[[i]] <- multi_int_fun(i, 365*3, parms, init, amr, agg_func, ode_wrapper)
}

# Obtain the Integrals ----------------------------------------------------

opt_crit = data.frame("Flat_All" = integral(testrun_flat[[1]][[1]], 3000, 0.5), 
                     "Single_1" = integral(single_list[[1]][[1]], 3000, 0.5),
                     "Single_2" = integral(single_list[[2]][[1]], 3000, 0.5),
                     "Single_3" = integral(single_list[[3]][[1]], 3000, 0.5),
                     "Diff_1" = integral(diff_tax_list[[1]][[1]], 3000, 0.5),
                     "Diff_2" = integral(diff_tax_list[[2]][[1]], 3000, 0.5),
                     "Diff_3" = integral(diff_tax_list[[3]][[1]], 3000, 0.5),
                     "Diff_4" = integral(diff_tax_list[[4]][[1]], 3000, 0.5),
                     "Diff_5" = integral(diff_tax_list[[5]][[1]], 3000, 0.5),
                     "Diff_6" = integral(diff_tax_list[[6]][[1]], 3000, 0.5))

rownames(opt_crit) <- c("Tot_Inf", "AvgRes", "Avg_Anti", "Shannon")

rescale_data <- data.frame(t(apply(opt_crit[1:3,], MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
rescale_data["Shannon",] <- unlist(c("Flat_All" = NA, (opt_crit[4,2:10] - min(opt_crit[4,2:10]))/diff(range(opt_crit[4,2:10]))))

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

# Running the Baseline Model ----------------------------------------------

base <- ode_wrapper(y = init, func = amr, times = seq(0, 10000), parms = parms)[[1]]
base_totinf <- integral(base, 3000, 0.5)[1]
base_avgres_int <- integral(base, 3000, 0.5)[2]

# Identifying the Parameter Set -------------------------------------------

over_parms <- list(testrun_flat[[1]][[2]],
                   single_list[[1]][[2]],
                   single_list[[2]][[2]],
                   single_list[[3]][[2]],
                   diff_tax_list[[1]][[2]],
                   diff_tax_list[[2]][[2]],
                   diff_tax_list[[3]][[2]],
                   diff_tax_list[[4]][[2]],
                   diff_tax_list[[5]][[2]],
                   diff_tax_list[[6]][[2]])

list_usage <- list()

for(i in 1:10) {
  list_usage[[i]] <- usage_fun(over_parms[[i]])
}

# Relative Reduction Integral Comparison --------------------------------------------

reduc_usage_vec <- sapply(1:length(list_usage), function(x) sum(list_usage[[x]][,6]))
rel_AvgRes <- base_avgres_int - opt_crit[2,]
rel_totinf <- opt_crit[1,] - base_totinf

comp_totinf <- rel_totinf/reduc_usage_vec
comp_res <- rel_AvgRes/reduc_usage_vec

data.frame.rel <- rbind(comp_totinf, comp_res); rownames(data.frame.rel) <- c("Total Infections", "Average Resistance")
data.frame.rel["Available Antibiotics",] <- opt_crit["Avg_Anti",]
data.frame.rel["Shannon Index",] <- opt_crit["Shannon",]

rescale_data_rel <- data.frame(t(apply(data.frame.rel, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
rescale_data_rel["Available Antibiotics",] <- rescale_data["Avg_Anti",]
rescale_data_rel["Shannon Index",] <- rescale_data["Shannon",]

OG_melt_rel <- melt(as.matrix(data.frame.rel), measure.vars = colnames(data.frame.rel))
OG_melt_rel$Var2 <- factor(rep(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)", "Single Tax (LR)",
                               "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)", "Diff Tax (4 Rd)",
                               "Diff Tax (5 Rd)","Diff Tax (6 Rd)"), each = 4))

OG_melt_rel$Var2  <- factor(OG_melt_rel$Var2 , levels = unique(OG_melt_rel$Var2))

rescale_melt_rel <- melt(as.matrix(rescale_data_rel), measure.vars = colnames(rescale_data_rel))

rescale_melt_rel$Var2 <- factor(rep(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)", "Single Tax (LR)",
                               "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)", "Diff Tax (4 Rd)",
                               "Diff Tax (5 Rd)","Diff Tax (6 Rd)"), each = 4))

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


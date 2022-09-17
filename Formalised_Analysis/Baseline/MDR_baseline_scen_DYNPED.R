library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
library("rbenchmark")

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


# Function to Remove NAs and Round Data -----------------------------------

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

# Function to Extract Antibiotic Usage Reductions -------------------------

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


# Flat Tax Model ----------------------------------------------------------

single_tax <- function(res_order, tax, parms, init, func, agg_func) {
  
  #First Run
  parms[["base_tax"]] <- tax
  
  run_1rd <- agg_func(remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))))
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  parms[["eff_tax"]][as.numeric(substr(res_order_vec, 2, 2)), c(1:6)] <- as.numeric(parms[["base_tax"]])
  
  #Real Model Run 
  run_real <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms)))
  return(run_real)
}


# Single Taxation Function ------------------------------------------------

single_tax <- function(res_order, tax, parms, init, func, agg_func) {
  
  #First Run
  parms[["base_tax"]] <- tax
  
  run_1rd <- agg_func(remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))))
  values_1rd <- tail(run_1rd, 1)[4:6]
  
  res_order_vec <- c(names(values_1rd)[which.max(values_1rd)],
                     names(values_1rd)[setdiff(1:3, c(which.min(values_1rd), which.max(values_1rd)))],
                     names(values_1rd)[which.min(values_1rd)])[res_order]
  
  parms[["eff_tax"]][as.numeric(substr(res_order_vec, 2, 2)), c(1:6)] <- as.numeric(parms[["base_tax"]])
  
  #Real Model Run 
  run_real <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms)))
  return(run_real)
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
             PED = matrix(c(-1, 0, 0, 
                            0, -1, 0,
                            0, 0, -1), #Be aware of this matrix
                          nrow = 3, ncol = 3, byrow = T),
             eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0, 
                                0, 0, 0, 0, 0, 0), 
                              nrow = 3, ncol = 6, byrow = T),
             t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5)

# Testing the Outcome Measure ---------------------------------------------

parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5

#Baseline Test Run
testrun_flat <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
integral(testrun_flat, 3000, 0.75) 

#Baseline Differential Taxation Run 
diff_test_run <- multi_int_fun(6, 365*3, parms, init, amr, agg_func)[[1]]
integral(diff_test_run, 3000, 0.75) 

# Baseline Model ----------------------------------------------------------

parms1 <- parms
testrun_flat <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
test_run_agg <- agg_func(testrun_flat) 

test_run_agg$AverageRes <- rowMeans(testrun_flat[,4:6])
test_run_agg$TotInf <- rowSums(testrun_flat[,3:10])

test_plot_flat <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[4:6])
ggplot(test_plot_flat, aes(time, value, color = variable)) + geom_line() + theme_bw()

# Intervention Scenarios --------------------------------------------------

#Flat Tax

parms1 <- parms; parms1[["eff_tax"]][,] <- 0.5
testrun_flat <- list(remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))))

#Single Tax 
single_list <- list()

for(i in 1:3) {
  parms1 <- parms
  single_list[[i]] <- data.frame(single_tax(i, 0.5, parms1, init, amr, agg_func),
                                 "scen" = paste0("single_", "eff_tax", i))
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  dat <- multi_int_fun(i, 365*3, parms, init, amr, agg_func)[[1]]
  diff_tax_list[[i]] <- data.frame(dat,
                                   "scen" = paste0("diff", i))
}

# Plotting the Scenarios --------------------------------------------------

#Create a combined list of all the scenarios
list_scen <- unlist(list(testrun_flat, single_list, diff_tax_list), recursive = FALSE)

melt_data <- list()

#Melt each one
for(i in 1:length(list_scen)) {
  data_agg <- agg_func(list_scen[[i]]) 
  colnames(data_agg)[4:6] <- c("High Res (HR)", "Medium Res (MR)", "Low Res (LR)") 
  melt_data[[i]] <- data.frame(melt(data_agg, id.vars = "time", measure.vars = colnames(data_agg)[4:6]),
                               "scen" = c("flat", "single1", "single2", "single3",
                                          "diff1", "diff2", "diff3", "diff4", "diff5", "diff6")[i])
}

p_data <- list()

#Plotting Loop
for(i in 1:length(melt_data)) {
  data <- melt_data[[i]]
  p_data[[i]] <- ggplot(data, aes(time, value, color = variable)) + geom_line() + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "Time", y = "Prevalence", title = c("Flat Tax",
                                                 "Single Tax (HR)","Single Tax (MR)","Single Tax (LR)",
                                                 "Diff Tax (1 Rd)", "Diff Tax (2 Rd)", "Diff Tax (3 Rd)",
                                                 "Diff Tax (4 Rd)", "Diff Tax (5 Rd)", "Diff Tax (6 Rd)")[i], color = "")
}

ggarrange(p_data[[1]], "", "",
          p_data[[2]], p_data[[3]],p_data[[4]], 
          p_data[[5]], p_data[[6]], p_data[[7]],
          p_data[[8]], p_data[[9]], p_data[[10]],
          labels = c("A", "", "",
                     "B", "", "",
                     "C", "", "",
                     "", "", ""), hjust = -.1,  nrow = 4, ncol = 3, common.legend = T, legend = "bottom")

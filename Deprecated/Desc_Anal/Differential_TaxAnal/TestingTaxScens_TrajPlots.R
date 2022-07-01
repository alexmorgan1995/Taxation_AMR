library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:6]),
                      sum(data_temp[4:6]),
                      sum(data_temp[3]),
                      sum(rowMeans(data_temp[4:6])),
                      sum(data_temp[4]),
                      sum(data_temp[5]),
                      sum(data_temp[6])),5)
  return(out_vec)
}

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    sigma_base1 <- sigma_1
    sigma_base2 <- sigma_2
    sigma_base3 <- sigma_3

    if(t > t_n) {
      sigma_1 <- sigma_base1*eff_tax1_1 
      sigma_2 <- sigma_base2*eff_tax2_1
      sigma_3 <- sigma_base3*eff_tax3_1
    }
    
    if(t > (t_n + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_2 
      sigma_2 <- sigma_base2*eff_tax2_2
      sigma_3 <- sigma_base3*eff_tax3_2
    }
    
    if(t > (t_n + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_3 
      sigma_2 <- sigma_base2*eff_tax2_3
      sigma_3 <- sigma_base3*eff_tax3_3
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_4 
      sigma_2 <- sigma_base2*eff_tax2_4
      sigma_3 <- sigma_base3*eff_tax3_4
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_5 
      sigma_2 <- sigma_base2*eff_tax2_5
      sigma_3 <- sigma_base3*eff_tax3_5
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_6 
      sigma_2 <- sigma_base2*eff_tax2_6
      sigma_3 <- sigma_base3*eff_tax3_6
    }

        dX = - beta*X*(W+(R1*c1)+(R2*c2)+(R3*c3)) + (1-sigma_1+sigma_2+sigma_3)*mu_w*W + (sigma_1+sigma_2+sigma_3)*mu_t*W +
      (1-sigma_2+sigma_3)*mu_r*R1 + (sigma_2+sigma_3)*mu_t*R1 +
      (1-sigma_1+sigma_3)*mu_r*R2 + (sigma_1+sigma_3)*mu_t*R2 + 
      (1-sigma_1+sigma_2)*mu_r*R3 + (sigma_1+sigma_2)*mu_t*R3
    
    dW = beta*X*W - (1-sigma_1+sigma_2+sigma_3)*mu_w*W - (sigma_1+sigma_2+sigma_3)*mu_t*W + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c1*R1 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c2*R2 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c3*R3 - 
      sigma_1*eta_wr*W - sigma_2*eta_wr*W - sigma_3*eta_wr*W
    
    dR1 = beta*X*R1*c1 - (1-sigma_2+sigma_3)*mu_r*R1 - (sigma_2+sigma_3)*mu_t*R1 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c1*R1 + sigma_1*eta_wr*W
    
    dR2 = beta*X*R2*c2 - (1-sigma_1+sigma_3)*mu_r*R2 - (sigma_1+sigma_3)*mu_t*R2 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c2*R2 + sigma_2*eta_wr*W
    
    dR3 = beta*X*R3*c3 - (1-sigma_1+sigma_2)*mu_r*R3 - (sigma_1+sigma_2)*mu_t*R3 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c3*R3 + sigma_3*eta_wr*W
    
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX,dW,dR1,dR2,dR3)))
  })
}


# Dual Model --------------------------------------------------------------

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 = 0, R3 = 0)

multi_int_fun <- function(int_gen, basetax, time_between, parms, init, func){

  parms["time_between"] <- time_between
    
    #First Run
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))
    values <- tail(testrun, 1)[4:6]
    
    parms[grep("eff_tax1", names(parms), value = T)] <- 1-(basetax*(values[1]/values[2]))
    parms[grep("eff_tax2", names(parms), value = T)] <- 1-(basetax*(values[2]/values[2]))
    parms[grep("eff_tax3", names(parms), value = T)] <- 1-(basetax*(values[3]/values[2]))
    parms[parms < 0] <- 0

    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    
    #Second Run
    if(int_gen >= 2){
      testrun1 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms))
      values1 <- tail(testrun1, 1)[4:6]
      
      low_char <- names(values1)[which.min(values1)]
      high_char <- names(values1)[which.max(values1)]
      med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]

      parms[grep(paste0("eff_tax",substr(low_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[low_char],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[med_char],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[high_char],1)/values[2]))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 3){
      testrun2 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms))
      values2 <- tail(testrun2, 1)[4:6]
      
      low_char1 <- names(values2)[which.min(values2)]
      high_char1 <- names(values2)[which.max(values2)]
      med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
      
      parms[grep(paste0("eff_tax",substr(low_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[low_char1],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[med_char1],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[high_char1],1)/values[2]))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 4){
      testrun3 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms))
      values3 <- tail(testrun3, 1)[4:6]
      
      low_char2 <- names(values3)[which.min(values3)]
      high_char2 <- names(values3)[which.max(values3)]
      med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]

      parms[grep(paste0("eff_tax",substr(low_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[low_char2],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[med_char2],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[high_char2],1)/values[2]))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 5){
      testrun4 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms))
      values4 <- tail(testrun4, 1)[4:6]
      
      low_char3 <- names(values4)[which.min(values4)]
      high_char3 <- names(values4)[which.max(values4)]
      med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
      
      parms[grep(paste0("eff_tax",substr(low_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[low_char3],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[med_char3],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[high_char3],1)/values[2]))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 6){
      testrun5 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms))
      values5 <- tail(testrun5, 1)[4:6]
      
      low_char4 <- names(values5)[which.min(values5)]
      high_char4 <- names(values5)[which.max(values5)]
      med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
      
      parms[grep(paste0("eff_tax",substr(low_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[low_char4],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[med_char4],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[high_char4],1)/values[2]))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    
    return(list(testrun, parms))
}

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 = 0, R3 = 0)

parms = c(beta = 5, sigma_1 = 0.25, sigma_2 = 0.25, sigma_3 = 0.25,
          mu_w = 1/12, mu_r = 1/10, 
          mu_t = 1/7, eta_wr = 0.01, eta_rw= 0.01, 
          c1= 0.95, c2 = 0.94, c3 = 0.91,
          eff_tax1_1 = 1, eff_tax2_1 = 1, eff_tax3_1 = 1, 
          eff_tax1_2 = 1, eff_tax2_2 = 1, eff_tax3_2 = 1, 
          eff_tax1_3 = 1, eff_tax2_3 = 1, eff_tax3_3 = 1, 
          eff_tax1_4 = 1, eff_tax2_4 = 1, eff_tax3_4 = 1, 
          eff_tax1_5 = 1, eff_tax2_5 = 1, eff_tax3_5 = 1, 
          eff_tax1_6 = 1, eff_tax2_6 = 1, eff_tax3_6 = 1, 
          t_n = 3000, time_between = Inf)

# Run the Scenarios -------------------------------------------------------

#Flat Tax Model - ALL

parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))

test <- melt(testrun_flat, id.vars = c("time"), measure.vars = c(colnames(testrun_flat)[3:6]))
ggplot(test, aes(time, value, color = as.character(variable))) + geom_line()

#Dual Taxation 

dual_p <- matrix(c(1,1, 2,2, 3,3), ncol = 2, nrow = 3) 

dual_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", dual_p[i,][1], "|", "eff_tax", dual_p[i,][2]), names(parms1), value = TRUE)] <- 0.5
  dual_list[[i]] <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1),
                               scen = paste0("dualtax", dual_p[i,][1], "_", dual_p[i,][2]))
}

test <- melt(dual_list[[1]], id.vars = c("time"), measure.vars = c(colnames(testrun_flat)[3:6]))
ggplot(test, aes(time, value, color = as.character(variable))) + geom_line()

#Single Taxation 

single_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", i), names(parms1), value = TRUE)] <- 0.5
  single_list[[i]] <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1),
                                 "scen" = paste0("single_", "eff_tax", i))
}

test <- melt(single_list[[3]], id.vars = c("time"), measure.vars = c(colnames(testrun_flat)[3:6]))
ggplot(test, aes(time, value, color = as.character(variable))) + geom_line()

#Diff Taxation

diff_tax_list <- list()
for(i in 1:6) {
  diff_tax_list[[i]] <- data.frame(multi_int_fun(i, 0.5, 365*3, parms, init, amr)[[1]], 
                                   "scen" = paste0("diff", i))
}

test <- melt(diff_tax_list[[6]], id.vars = c("time"), measure.vars = c(colnames(testrun_flat)[3:6]))
ggplot(test, aes(time, value, color = as.character(variable))) + geom_line()

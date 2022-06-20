library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:6]),
                      sum(rowMeans(data_temp[4:6]))),5)
  return(out_vec)
}


# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    R_base1 <- R1
    R_base2 <- R2
    R_base3 <- R3
    
    if(t > t_n) {
      R1 <- R_base1*(1-eff_tax1_1) 
      R2 <- R_base2*(1-eff_tax2_1)
      R3 <- R_base3*(1-eff_tax3_1)
    }
    
    if(t > (t_n + time_between)) {
      R1 <- R_base1*(1-eff_tax1_2) 
      R2 <- R_base2*(1-eff_tax2_2)
      R3 <- R_base3*(1-eff_tax3_2)
    }
    
    if(t > (t_n + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_3) 
      R2 <- R_base2*(1-eff_tax2_3)
      R3 <- R_base3*(1-eff_tax3_3)
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_4) 
      R2 <- R_base2*(1-eff_tax2_4)
      R3 <- R_base3*(1-eff_tax3_4)
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_5) 
      R2 <- R_base2*(1-eff_tax2_5)
      R3 <- R_base3*(1-eff_tax3_5)
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      R1 <- R_base1*(1-eff_tax1_6) 
      R2 <- R_base2*(1-eff_tax2_6)
      R3 <- R_base3*(1-eff_tax3_6)
    }
    
  dS = lambda - lambda*S + 
    mu_wt*eta*Wc + mu_r1*eta*Rc1 + mu_r2*eta*Rc2 + mu_r3*eta*Rc3 + 
    P1*theta + P2*theta + P3*theta + 
    Pr*rho*(R1+R2+R3)*Wc + 
    Pr*(R2+R3)*Rc1 + Pr*(R1+R3)*Rc2 + Pr*(R1+R2)*Rc3 - 
    Pr*(R1+R2+R3)*S - 
    beta*S*(Wc + Wi) - beta*S*(Ri1 + Rc1)*fc1 - beta*S*(Ri2 + Rc2)*fc2 - beta*S*(Ri3 + Rc3)*fc3 +
    tau*(R1+R2+R3)*Wi + tau*(R2+R3)*Ri1 + tau*(R1+R3)*Ri2 + tau*(R1+R2)*Ri3   
    
  dWc = -lambda*Wc - mu_wt*eta*Wc + mu_wt*Wi - Pr*(R1+R2+R3)*Wc + 
    beta*S*(Wc + Wi)*psi + kappa*(Rc1 + Rc2 + Rc3) 
  
  dWi = -lambda*Wi- mu_wt*Wi - tau*(R1+R2+R3)*Wi + beta*S*(Wc + Wi)*(1-psi)
  
  dRc1 = -lambda*Rc1 - mu_r1*eta*Rc1 + mu_r1*Ri1 + Pr*R1*(1-rho)*Wc - Pr*(R2+R3)*Rc1 + 
    beta*S*(Ri1 + Rc1)*psi*fc1 - kappa*Rc1 + Pr*R1*(1-rho)*S
  
  dRi1 = -lambda*Ri1 - mu_r1*Ri1 - tau*(R2 + R3)*Ri1 + beta*S*(Ri1 + Rc1)*(1-psi)*fc1
  
  dRc2 = -lambda*Rc2 - mu_r2*eta*Rc2 + mu_r2*Ri2 + Pr*R2*(1-rho)*Wc - Pr*(R1+R3)*Rc2 + 
    beta*S*(Ri2 + Rc2)*psi*fc2 - kappa*Rc2 + Pr*R2*(1-rho)*S
  
  dRi2 = -lambda*Ri2 - mu_r2*Ri2 - tau*(R1 + R3)*Ri2 + beta*S*(Ri2 + Rc2)*(1-psi)*fc2
  
  dRc3 = -lambda*Rc3 - mu_r3*eta*Rc3 + mu_r3*Ri3 + Pr*R3*(1-rho)*Wc - Pr*(R1+R2)*Rc3 + 
    beta*S*(Ri3 + Rc3)*psi*fc3 - kappa*Rc3 + Pr*R3*(1-rho)*S
  
  dRi3 = -lambda*Ri3 - mu_r3*Ri3 - tau*(R1 + R2)*Ri3 + beta*S*(Ri3 + Rc3)*(1-psi)*fc3
  
  dP1 = -lambda*P1 + Pr*R1*rho*S - P1*theta
  
  dP2 = -lambda*P2 + Pr*R2*rho*S - P2*theta
  
  dP3 = -lambda*P3 + Pr*R3*rho*S - P3*theta
  
    #Calculating the Proportion Integrals
    
    return(list(c(dS, dWc, dWi, dRc1, dRi1, dRc2, dRi2, dRc3, dRi3, dP1, dP2, dP3)))
  })
}

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, basetax, time_between, parms, init, func){
  
  parms["time_between"] <- time_between
  
  #First Run
  testrun <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))
  values <- tail(testrun, 1)[4:6]
  
  parms[grep("eff_tax1", names(parms), value = T)] <- (basetax*(values[1]/values[2]))
  parms[grep("eff_tax2", names(parms), value = T)] <- (basetax*(values[2]/values[2]))
  parms[grep("eff_tax3", names(parms), value = T)] <- (basetax*(values[3]/values[2]))
  parms[parms < 0] <- 0
  
  testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  
  #Second Run
  if(int_gen >= 2){
    testrun1 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms))
    values1 <- tail(testrun1, 1)[4:6]
    
    low_char <- names(values1)[which.min(values1)]
    high_char <- names(values1)[which.max(values1)]
    med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]
    
    parms[grep(paste0("eff_tax",substr(low_char, 2, 2)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[low_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char, 2, 2)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[med_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char, 2, 2)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[high_char],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 3){
    testrun2 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms))
    values2 <- tail(testrun2, 1)[4:6]
    
    low_char1 <- names(values2)[which.min(values2)]
    high_char1 <- names(values2)[which.max(values2)]
    med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
    
    parms[grep(paste0("eff_tax",substr(low_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[low_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[med_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[high_char1],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 4){
    testrun3 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms))
    values3 <- tail(testrun3, 1)[4:6]
    
    low_char2 <- names(values3)[which.min(values3)]
    high_char2 <- names(values3)[which.max(values3)]
    med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]
    
    parms[grep(paste0("eff_tax",substr(low_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[low_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[med_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[high_char2],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 5){
    testrun4 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms))
    values4 <- tail(testrun4, 1)[4:6]
    
    low_char3 <- names(values4)[which.min(values4)]
    high_char3 <- names(values4)[which.max(values4)]
    med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
    
    parms[grep(paste0("eff_tax",substr(low_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[low_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[med_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[high_char3],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 6){
    testrun5 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms))
    values5 <- tail(testrun5, 1)[4:6]
    
    low_char4 <- names(values5)[which.min(values5)]
    high_char4 <- names(values5)[which.max(values5)]
    med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
    
    parms[grep(paste0("eff_tax",substr(low_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[low_char4],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[med_char4],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[high_char4],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  
  return(list(testrun, parms))
}

# Run the ODE -------------------------------------------------------------

#Init Parameters 
init <- c(S = 0.99, Wc = 1-0.99, Wi = 0, 
          Rc1 = 0, Ri1 = 0,
          Rc2 = 0, Ri2 = 0,
          Rc3 = 0, Ri3 = 0,
          P1 = 0, P2 = 0, P3 = 0)

parms = c(beta = 5, R1 = 0.25, R2 = 0.25, R3 = 0.25,
          lambda = 0.1, Pr = 1, rho = 0.8, theta = 0.1,
          tau = 1, psi = 0.95, kappa = 0.01, eta = 0.1, 
          fc1 = 0.9, fc2 = 0.86, fc3 = 0.7, 
          mu_wt = 0.01, mu_r1 = 0.01, mu_r2 = 0.01, mu_r3 = 0.01,  
          
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          t_n = 10, time_between = Inf)

#Flat Tax
parms1 <- parms
#parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 1000), parms = parms1))
#testrun_flat$AverageRes <- rowMeans(testrun_flat[,4:6])
#testrun_flat$TotInf <- rowSums(testrun_flat[,3:6])

testrun_flat$WT <- rowSums(testrun_flat[,3:4])
testrun_flat$Res1 <- rowSums(testrun_flat[,5:6])
testrun_flat$Res2 <- rowSums(testrun_flat[,7:8])
testrun_flat$Res3 <- rowSums(testrun_flat[,9:10])
testrun_flat$proph <- rowSums(testrun_flat[,11:13])

t_melt <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[-1])
t_melt <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[14:18])

ggplot(t_melt, aes(time, value, color = variable)) + geom_line()


# Run the Models ----------------------------------------------------------

#Flat Tax
parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))
testrun_flat$AverageRes <- rowMeans(testrun_flat[,4:6])
testrun_flat$TotInf <- rowSums(testrun_flat[,3:6])

#Dual Tax
dual_p <- matrix(c(1,1, 2,2, 3,3), ncol = 2, nrow = 3) 
dual_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", dual_p[i,][1], "|", "eff_tax", dual_p[i,][2]), names(parms1), value = TRUE)] <- 0.5
  dat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))
  dual_list[[i]] <- data.frame(dat,
                               "AverageRes" = rowMeans(dat[,4:6]),
                               "TotInf" = rowSums(dat[,3:6]),
                               scen = paste0("dualtax", dual_p[i,][1], "_", dual_p[i,][2]))
}

#Single Tax 
single_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", i), names(parms1), value = TRUE)] <- 0.5
  single_list[[i]] <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1),
                                 "scen" = paste0("single_", "eff_tax", i))
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  dat <- multi_int_fun(i, 0.5, 365*3, parms, init, amr)[[1]]
  diff_tax_list[[i]] <- data.frame(dat,
                                   "AverageRes" = rowMeans(dat[,4:6]),
                                   "TotInf" = rowSums(dat[,3:6]),
                                   "scen" = paste0("diff", i))
}

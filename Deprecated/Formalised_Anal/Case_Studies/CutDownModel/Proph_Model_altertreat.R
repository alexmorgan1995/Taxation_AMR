library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  data_temp <- data_temp[,-14]
  data_temp$Res1 <- rowSums(data_temp[,5:6])
  data_temp$Res2 <- rowSums(data_temp[,7:8])
  data_temp$Res3 <- rowSums(data_temp[,9:10])

  out_vec <- signif(c(sum(data_temp[3:10]),
                      sum(rowMeans(data_temp[14:16]))),
                    5)
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
    
    R1 <- ifelse(R1 > 0, R1, 0)
    R2 <- ifelse(R2 > 0, R2, 0)
    R3 <- ifelse(R3 > 0, R3, 0)
    
  dS = lambda - lambda*S + 
    mu_wt*eta*Wc + mu_r*eta*Rc1 + mu_r*eta*Rc2 + mu_r*eta*Rc3 + 
    P1*theta + P2*theta + P3*theta + 
    Pr*zeta*(R1+R2+R3)*Wc + 
    Pr*zeta*(R2+R3)*Rc1 + Pr*zeta*(R1+R3)*Rc2 + Pr*zeta*(R1+R2)*Rc3 - 
    Pr*(R1+R2+R3)*S - 
    beta*S*(Wc + Wi) - beta*S*(Ri1 + Rc1)*fc1 - beta*S*(Ri2 + Rc2)*fc2 - beta*S*(Ri3 + Rc3)*fc3 +
    tau*(R1+R2+R3)*Wi + tau*(R2+R3)*Ri1 + tau*(R1+R3)*Ri2 + tau*(R1+R2)*Ri3   
    
  dWc = -lambda*Wc - mu_wt*eta*Wc + mu_wt*Wi - Pr*zeta*(R1+R2+R3)*Wc + 
    beta*S*(Wc + Wi)*psi + kappa*(Rc1 + Rc2 + Rc3) 
  
  dWi = -lambda*Wi- mu_wt*Wi - tau*(R1+R2+R3)*Wi + beta*S*(Wc + Wi)*(1-psi)
  
  dRc1 = -lambda*Rc1 - mu_r*eta*Rc1 + mu_r*Ri1 + Pr*R1*(1-zeta)*Wc - Pr*zeta*(R2+R3)*Rc1 + #The last zeta term here corresponds to a reduce effect of treatment in Resistant infections 
    beta*S*(Ri1 + Rc1)*psi*fc1 - kappa*Rc1 + Pr*R1*(1-rho)*S #Ideally we would model treatment failure here - but that would require multidrug resistance or replacement of strains - so ignored the 1-zeta here 
  
  dRi1 = -lambda*Ri1 - mu_r*Ri1 - tau*(R2 + R3)*Ri1 + beta*S*(Ri1 + Rc1)*(1-psi)*fc1
  
  dRc2 = -lambda*Rc2 - mu_r*eta*Rc2 + mu_r*Ri2 + Pr*R2*(1-zeta)*Wc - Pr*zeta*(R1+R3)*Rc2 + 
    beta*S*(Ri2 + Rc2)*psi*fc2 - kappa*Rc2 + Pr*R2*(1-rho)*S
  
  dRi2 = -lambda*Ri2 - mu_r*Ri2 - tau*(R1 + R3)*Ri2 + beta*S*(Ri2 + Rc2)*(1-psi)*fc2
  
  dRc3 = -lambda*Rc3 - mu_r*eta*Rc3 + mu_r*Ri3 + Pr*R3*(1-zeta)*Wc - Pr*zeta*(R1+R2)*Rc3 + 
    beta*S*(Ri3 + Rc3)*psi*fc3 - kappa*Rc3 + Pr*R3*(1-rho)*S
  
  dRi3 = -lambda*Ri3 - mu_r*Ri3 - tau*(R1 + R2)*Ri3 + beta*S*(Ri3 + Rc3)*(1-psi)*fc3
  
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
  testrun$Res1 <- rowSums(testrun[,5:6])
  testrun$Res2 <- rowSums(testrun[,7:8])
  testrun$Res3 <- rowSums(testrun[,9:10])
  values <- tail(testrun, 1)[14:16]
  
  parms[grep("eff_tax1", names(parms), value = T)] <- (basetax*(values[1]/values[2]))
  parms[grep("eff_tax2", names(parms), value = T)] <- (basetax*(values[2]/values[2]))
  parms[grep("eff_tax3", names(parms), value = T)] <- (basetax*(values[3]/values[2]))
  parms[parms < 0] <- 0
  
  testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  
  #Second Run
  if(int_gen >= 2){
    testrun1 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms))
    testrun1$Res1 <- rowSums(testrun1[,5:6])
    testrun1$Res2 <- rowSums(testrun1[,7:8])
    testrun1$Res3 <- rowSums(testrun1[,9:10])
    values1 <- tail(testrun1, 1)[14:16]
    
    low_char <- names(values1)[which.min(values1)]
    high_char <- names(values1)[which.max(values1)]
    med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]
    

    parms[grep(paste0("eff_tax",substr(low_char, 4, 4)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[low_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char, 4, 4)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[med_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char, 4, 4)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[high_char],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 3){
    testrun2 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms))
    testrun2$Res1 <- rowSums(testrun2[,5:6])
    testrun2$Res2 <- rowSums(testrun2[,7:8])
    testrun2$Res3 <- rowSums(testrun2[,9:10])
    values2 <- tail(testrun2, 1)[14:16]
    
  
    low_char1 <- names(values2)[which.min(values2)]
    high_char1 <- names(values2)[which.max(values2)]
    med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
    
    
    parms[grep(paste0("eff_tax",substr(low_char1, 4, 4)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[low_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char1, 4, 4)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[med_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char1, 4, 4)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[high_char1],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 4){
    testrun3 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms))
    testrun3$Res1 <- rowSums(testrun3[,5:6])
    testrun3$Res2 <- rowSums(testrun3[,7:8])
    testrun3$Res3 <- rowSums(testrun3[,9:10])

    values3 <- tail(testrun3, 1)[14:16]
    
    low_char2 <- names(values3)[which.min(values3)]
    high_char2 <- names(values3)[which.max(values3)]
    med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]
    
    parms[grep(paste0("eff_tax",substr(low_char2, 4, 4)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[low_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char2, 4, 4)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[med_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char2, 4, 4)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[high_char2],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 5){
    testrun4 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms))
    testrun4$Res1 <- rowSums(testrun4[,5:6])
    testrun4$Res2 <- rowSums(testrun4[,7:8])
    testrun4$Res3 <- rowSums(testrun4[,9:10])
    values4 <- tail(testrun4, 1)[14:16]
    
    low_char3 <- names(values4)[which.min(values4)]
    high_char3 <- names(values4)[which.max(values4)]
    med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
    
    parms[grep(paste0("eff_tax",substr(low_char3, 4, 4)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[low_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char3, 4, 4)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[med_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char3, 4, 4)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[high_char3],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
  }
  
  if(int_gen >= 6){
    testrun5 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms))
    testrun5$Res1 <- rowSums(testrun5[,5:6])
    testrun5$Res2 <- rowSums(testrun5[,7:8])
    testrun5$Res3 <- rowSums(testrun5[,9:10])
    values5 <- tail(testrun5, 1)[14:16]
    
    low_char4 <- names(values5)[which.min(values5)]
    high_char4 <- names(values5)[which.max(values5)]
    med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
    
    parms[grep(paste0("eff_tax",substr(low_char4, 4, 4)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[low_char4],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char4, 4, 4)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[med_char4],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char4, 4, 4)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[high_char4],1)/values[2]))
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
          mu_wt = 1/12, mu_r = 1/10,  
          
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          t_n = 1000, time_between = Inf)


#Flat Tax
parms1 <- parms
#parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))
#testrun_flat$AverageRes <- rowMeans(testrun_flat[,4:6])
#testrun_flat$TotInf <- rowSums(testrun_flat[,3:6])

testrun_flat$WT <- rowSums(testrun_flat[,3:4])
testrun_flat$Res1 <- rowSums(testrun_flat[,5:6])
testrun_flat$Res2 <- rowSums(testrun_flat[,7:8])
testrun_flat$Res3 <- rowSums(testrun_flat[,9:10])
testrun_flat$proph <- rowSums(testrun_flat[,11:13])

t_melt <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[14:18])
#t_melt <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[-1])

ggplot(t_melt, aes(time, value, color = variable)) + geom_line(size = 1.2) + theme_bw() + 
  scale_x_continuous(name = "Time (days)", expand = c(0, 0)) +  scale_y_continuous(name = "Proportion Infected", expand = c(0, 0), limits = c(0,0.5)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") 

# Run the Models ----------------------------------------------------------

parms1 = c(beta = 5, R1 = 0.25, R2 = 0.25, R3 = 0.25,
          lambda = 0.1, Pr = 1, rho = 0.8, theta = 0.1,
          tau = 1, psi = 0.95, kappa = 0.01, eta = 0.1, 
          fc1 = 0.9, fc2 = 0.86, fc3 = 0.7, 
          mu_wt = 1/12, mu_r = 1/10,  
          
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          t_n = 1000, time_between = Inf)


#Flat Tax
parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))

testrun_flat$Res1 <- rowSums(testrun_flat[,5:6])
testrun_flat$Res2 <- rowSums(testrun_flat[,7:8])
testrun_flat$Res3 <- rowSums(testrun_flat[,9:10])

t_melt <- melt(testrun_flat, id.vars = "time", measure.vars = c("Res1", "Res2", "Res3"))

ggplot(t_melt, aes(time, value, color = variable)) + geom_line()

#Dual Tax
dual_p <- matrix(c(1,1, 2,2, 3,3), ncol = 2, nrow = 3) 
dual_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", dual_p[i,][1], "|", "eff_tax", dual_p[i,][2]), names(parms1), value = TRUE)] <- 0.5
  dat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1))
  dual_list[[i]] <- data.frame(dat,
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
                                   "scen" = paste0("diff", i))
}



test <- multi_int_fun(5, 0.5, 365*3, parms, init, amr)[[1]]
test$WT <- rowSums(test[,3:4])
test$Res1 <- rowSums(test[,5:6])
test$Res2 <- rowSums(test[,7:8])
test$Res3 <- rowSums(test[,9:10])
test$proph <- rowSums(test[,11:13])

m_test <- melt(test, id.vars = "time", measure.vars = c("WT","Res1", "Res2", "Res3", "proph"))

ggplot(m_test, aes(time, value, color = variable)) + geom_line(size = 1.2) + theme_bw() + 
  scale_x_continuous(name = "Time (days)", expand = c(0, 0)) +  scale_y_continuous(name = "Proportion Infected", expand = c(0, 0), limits = c(0,1)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") 

# Obtain the Integrals ----------------------------------------------------

inf_int = data.frame("Flat_All" = integral(testrun_flat, 1000), 
                     "Dual_12" = integral(dual_list[[1]], 1000), 
                     "Dual_13" = integral(dual_list[[2]], 1000), 
                     "Dual_23" = integral(dual_list[[3]], 1000), 
                     "Single_1" = integral(single_list[[1]], 1000),
                     "Single_2" = integral(single_list[[2]], 1000),
                     "Single_3" = integral(single_list[[3]], 1000),
                     "Diff_1" = integral(diff_tax_list[[1]], 1000),
                     "Diff_2" = integral(diff_tax_list[[2]], 1000),
                     "Diff_3" = integral(diff_tax_list[[3]], 1000),
                     "Diff_4" = integral(diff_tax_list[[4]], 1000),
                     "Diff_5" = integral(diff_tax_list[[5]], 1000),
                     "Diff_6" = integral(diff_tax_list[[6]], 1000))

rownames(inf_int) <- c("Total Infections", "Average Resistance")


# Integrals ---------------------------------------------------------------

diff_tax_list <- list()

for(i in 1:6) {
  diff_tax_list[[i]] <- data.frame(multi_int_fun(i, 0.5, 365*3, parms, init, amr)[[2]], 
                                   "scen" = paste0("diff", i))
}

#Find The Resistance Lost 
#Baseline

base <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms))
base_totinf <- integral(base, 1000)[1]
base_avgres_int <- integral(base, 1000)[2]

#Lost Resistance Function 
#Supply a dataframe with time and the effective tax over time for class 1, 2 and 3.
#I need to extract the 1-6 data with the 

#Change in Usage From Baseline - Only for differential taxation models 

usage_fun <- function(parms){
  usage = data.frame("time" = seq(0,9000),
                     "PopUsage1" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax1_", x)]])*0.25, 365*3)),
                                     rep((1-parms[["eff_tax1_6"]])*0.25, 9001-365*3*6)),
                     "PopUsage2" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax2_", x)]])*0.25, 365*3)),
                                     rep((1-parms[["eff_tax2_6"]])*0.25, 9001-365*3*6)),
                     "PopUsage3" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax3_", x)]])*0.25, 365*3)),
                                     rep((1-parms[["eff_tax3_6"]])*0.25, 9001-365*3*6)))
  usage[usage < 0] <- 0
  usage$totusage = rowSums(usage[2:4])
  usage$reduc_use = 0.75 - usage$totusage
  return(usage)
}

#Now I need to do this for all 6 scenarios. 
parms_flat <- parms; parms_flat[grep("eff_tax", names(parms), value =T)]  <- 0.5
parms_single1 <- parms; parms_single1[grep("eff_tax1", names(parms), value =T)]  <- 0.5
parms_single2 <- parms; parms_single2[grep("eff_tax2", names(parms), value =T)]  <- 0.5
parms_single3 <- parms; parms_single3[grep("eff_tax3", names(parms), value =T)]  <- 0.5
parms_dual12 <- parms; parms_dual12[grep("eff_tax1|eff_tax2", names(parms), value =T)]  <- 0.5
parms_dual23 <- parms; parms_dual23[grep("eff_tax2|eff_tax3", names(parms), value =T)]  <- 0.5
parms_dual13 <- parms; parms_dual13[grep("eff_tax1|eff_tax2", names(parms), value =T)]  <- 0.5

test <- usage_fun(parms_dual12)
sum(test[,6])

parm_list <- list(parms_flat, parms_single1, parms_single2, parms_single3,
                  parms_dual12, parms_dual23, parms_dual13)

over_parms <- append(parm_list, diff_tax_list)
list_usage <- lapply(over_parms, usage_fun)

# Relative Reduction Integral Comparison --------------------------------------------

reduc_usage_vec <- sapply(1:length(list_usage), function(x) sum(list_usage[[x]][,6]))
rel_AvgRes <- base_avgres_int - inf_int[2,]
rel_totinf <- inf_int[1,] - base_totinf

comp_totinf <- rel_totinf/reduc_usage_vec
comp_res <- rel_AvgRes/reduc_usage_vec

data.frame.rel <- rbind(comp_totinf, comp_res); rownames(data.frame.rel) <- c("Total Infections", "Average Resistance")

rescale_data_rel <- t(apply(data.frame.rel, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
OG_melt_rel <- melt(as.matrix(data.frame.rel), measure.vars = colnames(data.frame.rel))
rescale_melt_rel <- melt(as.matrix(rescale_data_rel), measure.vars = colnames(rescale_data_rel))
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
        strip.text = element_blank(), legend.position="none")

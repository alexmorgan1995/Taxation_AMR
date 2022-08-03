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
      sigma_1 <- sigma_base1*(1-(eff_tax1_1*PED1)) 
      sigma_2 <- sigma_base2*(1-(eff_tax2_1*PED2))
      sigma_3 <- sigma_base3*(1-(eff_tax3_1*PED3))
    }
    
    if(t > (t_n + time_between)) {
      sigma_1 <- sigma_base1*(1-(eff_tax1_2*PED1)) 
      sigma_2 <- sigma_base2*(1-(eff_tax2_2*PED2)) 
      sigma_3 <- sigma_base3*(1-(eff_tax3_2*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between)) {
      sigma_1 <- sigma_base1*(1-(eff_tax1_3*PED1)) 
      sigma_2 <- sigma_base2*(1-(eff_tax2_3*PED2)) 
      sigma_3 <- sigma_base3*(1-(eff_tax3_3*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*(1-(eff_tax1_4*PED1)) 
      sigma_2 <- sigma_base2*(1-(eff_tax2_4*PED2)) 
      sigma_3 <- sigma_base3*(1-(eff_tax3_4*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*(1-(eff_tax1_5*PED1)) 
      sigma_2 <- sigma_base2*(1-(eff_tax2_5*PED2)) 
      sigma_3 <- sigma_base3*(1-(eff_tax3_5*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*(1-(eff_tax1_6*PED1)) 
      sigma_2 <- sigma_base2*(1-(eff_tax2_6*PED2)) 
      sigma_3 <- sigma_base3*(1-(eff_tax3_6*PED3)) 
    }
    
    sigma_1 <- ifelse(sigma_1 > 0, sigma_1, 0)
    sigma_2 <- ifelse(sigma_2 > 0, sigma_2, 0)
    sigma_3 <- ifelse(sigma_3 > 0, sigma_3, 0)
    
    dX = - beta*X*(W+(R1*c1)+(R2*c2)+(R3*c3)) + (1-sigma_1+sigma_2+sigma_3)*mu_w*W + (sigma_1+sigma_2+sigma_3)*mu_t*W*(1-rho) +
      (1-sigma_2+sigma_3)*mu_r*R1 + (sigma_2+sigma_3)*mu_t*R1 +
      (1-sigma_1+sigma_3)*mu_r*R2 + (sigma_1+sigma_3)*mu_t*R2 + 
      (1-sigma_1+sigma_2)*mu_r*R3 + (sigma_1+sigma_2)*mu_t*R3
    
    dW = beta*X*W - (1-sigma_1+sigma_2+sigma_3)*mu_w*W - (sigma_1+sigma_2+sigma_3)*mu_t*W*(1-rho) + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R1 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R2 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R3 - 
      sigma_1*eta_wr*W*rho - sigma_2*eta_wr*W*rho - sigma_3*eta_wr*W*rho
    
    dR1 = beta*X*R1*c1 - (1-sigma_2+sigma_3)*mu_r*R1 - (sigma_2+sigma_3)*mu_t*R1 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R1 + sigma_1*eta_wr*W*rho
    
    dR2 = beta*X*R2*c2 - (1-sigma_1+sigma_3)*mu_r*R2 - (sigma_1+sigma_3)*mu_t*R2 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R2 + sigma_2*eta_wr*W*rho
    
    dR3 = beta*X*R3*c3 - (1-sigma_1+sigma_2)*mu_r*R3 - (sigma_1+sigma_2)*mu_t*R3 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*R3 + sigma_3*eta_wr*W*rho
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX,dW,dR1,dR2,dR3)))
  })
}


# Clean Model Run --------------------------------------------------------------

remNA_func <- function(dataframe){
  timing <- t(sapply(1:5, function(x)  dataframe[max(which(!is.na(dataframe[,x+1]))),]))
  if(timing[1,1] != 10000) {
    dataframe[seq(timing[[1]]+2,10001),2] <- timing[1,2]
    dataframe[seq(timing[[1]]+2,10001),3] <- timing[2,3]
    dataframe[seq(timing[[1]]+2,10001),4] <- timing[3,4]
    dataframe[seq(timing[[1]]+2,10001),5] <- timing[4,5]
    dataframe[seq(timing[[1]]+2,10001),6] <- timing[5,6]
    dataframe[dataframe < 1e-10] <- 0
  }
  return(dataframe)
}

# Dual Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, time_between, parms, init, func){
  
  parms["time_between"] <- time_between
  
  #First Run
  testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms, hmax = 1)))
  values <- tail(testrun, 1)[4:6]
  
  parms[grep("eff_tax1", names(parms), value = T)] <- (parms[["base_tax"]]*(values[1]/values[2]))
  parms[grep("eff_tax2", names(parms), value = T)] <- (parms[["base_tax"]]*(values[2]/values[2]))
  parms[grep("eff_tax3", names(parms), value = T)] <- (parms[["base_tax"]]*(values[3]/values[2]))
  parms[parms < 0] <- 0
  
  if(int_gen <= 2 | int_gen >= 6) {
    testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
  }
  
  #Second Run
  if(int_gen >= 2){
    
    testrun1 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms, hmax = 1)))
    values1 <- tail(testrun1, 1)[4:6]
    
    if(values1[1] == 0 & values1[2] == 0 & values1[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,3)]] <- 0
    } else {
      low_char <- names(values1)[which.min(values1)]
      high_char <- names(values1)[which.max(values1)]
      med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]
      parms[grep(paste0("eff_tax", substr(low_char, 2, 2)), names(parms), value = T)[-c(1)]] <- (parms[["base_tax"]]*(tail(testrun1[low_char],1)/values[2]))
      parms[grep(paste0("eff_tax", substr(med_char, 2, 2)), names(parms), value = T)[-c(1)]] <- (parms[["base_tax"]]*(tail(testrun1[med_char],1)/values[2]))
      parms[grep(paste0("eff_tax", substr(high_char, 2, 2)), names(parms), value = T)[-c(1)]] <- (parms[["base_tax"]]*(tail(testrun1[high_char],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 2) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 3){
    testrun2 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms, hmax = 1)))
    values2 <- tail(testrun2, 1)[4:6]
    
    if(values2[1] == 0 & values2[2] == 0 & values2[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,6)]] <- 0
    } else {
      low_char1 <- names(values2)[which.min(values2)]
      high_char1 <- names(values2)[which.max(values2)]
      med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
      parms[grep(paste0("eff_tax",substr(low_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (parms[["base_tax"]]*(tail(testrun2[low_char1],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (parms[["base_tax"]]*(tail(testrun2[med_char1],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- (parms[["base_tax"]]*(tail(testrun2[high_char1],1)/values[2]))
      
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 3) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 4){
    testrun3 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms, hmax = 1)))
    values3 <- tail(testrun3, 1)[4:6]
    
    if(values3[1] == 0 & values3[2] == 0 & values3[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,9)]] <- 0
    } else {
      low_char2 <- names(values3)[which.min(values3)]
      high_char2 <- names(values3)[which.max(values3)]
      med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]
      parms[grep(paste0("eff_tax",substr(low_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (parms[["base_tax"]]*(tail(testrun3[low_char2],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (parms[["base_tax"]]*(tail(testrun3[med_char2],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- (parms[["base_tax"]]*(tail(testrun3[high_char2],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 4) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 5){
    testrun4 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms, hmax = 1)))
    
    values4 <- tail(testrun4, 1)[4:6]
    
    if(values4[1] == 0 & values4[2] == 0 & values4[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,12)]] <- 0
    } else {
      low_char3 <- names(values4)[which.min(values4)]
      high_char3 <- names(values4)[which.max(values4)]
      med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
      parms[grep(paste0("eff_tax",substr(low_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (parms[["base_tax"]]*(tail(testrun4[low_char3],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (parms[["base_tax"]]*(tail(testrun4[med_char3],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- (parms[["base_tax"]]*(tail(testrun4[high_char3],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 5) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  if(int_gen >= 6){
    testrun5 <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms, hmax = 1)))
    values5 <- tail(testrun5, 1)[4:6]
    
    if(values5[1] == 0 & values5[2] == 0 & values5[3] == 0) {
      parms[grep(paste0("eff_tax"), names(parms), value = T)[-seq(1,15)]] <- 0
    } else {
      low_char4 <- names(values5)[which.min(values5)]
      high_char4 <- names(values5)[which.max(values5)]
      med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
      parms[grep(paste0("eff_tax",substr(low_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (parms[["base_tax"]]*(tail(testrun5[low_char4],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(med_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (parms[["base_tax"]]*(tail(testrun5[med_char4],1)/values[2]))
      parms[grep(paste0("eff_tax",substr(high_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (parms[["base_tax"]]*(tail(testrun5[high_char4],1)/values[2]))
      parms[parms < 0] <- 0
    }
    
    if(int_gen == 6) {
      testrun <- remNA_func(data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms, hmax = 1)))
    }
  }
  
  return(list(testrun, parms))
}

# Extract Usage -----------------------------------------------------------

usage_fun <- function(parms){
  usage = data.frame("time" = seq(0,7000),
                     "PopUsage1" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax1_", x)]])*parms[["sigma_1"]], 365*3)),
                                     rep((1-parms[["eff_tax1_6"]])*parms[["sigma_1"]], 7001-365*3*6)),
                     "PopUsage2" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax2_", x)]])*parms[["sigma_2"]], 365*3)),
                                     rep((1-parms[["eff_tax2_6"]])*parms[["sigma_2"]], 7001-365*3*6)),
                     "PopUsage3" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax3_", x)]])*parms[["sigma_3"]], 365*3)),
                                     rep((1-parms[["eff_tax3_6"]])*parms[["sigma_3"]], 7001-365*3*6)))
  usage[usage < 0] <- 0
  usage$totusage = rowSums(usage[2:4])
  usage$reduc_use = parms[["sigma_1"]] + parms[["sigma_2"]] + parms[["sigma_3"]] - usage$totusage
  return(usage)
}

# Baseline Parms ----------------------------------------------------------

init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 = 0, R3 = 0)

parms = c(beta = 5, sigma_1 = 0.25, sigma_2 = 0.25, sigma_3 = 0.25,
          mu_w = 1/12, mu_r = 1/10, 
          mu_t = 1/7, eta_wr = 0.15, eta_rw = 0.01, 
          c1 = 0.95, c2 = 0.91, c3 = 0.8,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 3000, time_between = Inf, rho = 0.5, base_tax = 0.5)

# Run the Models ----------------------------------------------------------

#Flat Tax
parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
testrun_flat$AverageRes <- rowMeans(testrun_flat[,4:6])
testrun_flat$TotInf <- rowSums(testrun_flat[,3:6])

test_plot_flat <- melt(testrun_flat, id.vars = "time", measure.vars = colnames(testrun_flat)[4:6])

ggplot(test_plot_flat, aes(time, value, color = variable)) + geom_line() + theme_bw()

#Dual Tax
dual_p <- matrix(c(1,1, 2,2, 3,3), ncol = 2, nrow = 3) 
dual_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", dual_p[i,][1], "|", "eff_tax", dual_p[i,][2]), names(parms1), value = TRUE)] <- 0.5
  dat <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)))
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
  single_list[[i]] <- data.frame(remNA_func(ode(y = init, func = amr, times = seq(0, 10000), parms = parms1)),
                                 "scen" = paste0("single_", "eff_tax", i))
}

#Diff Tax
diff_tax_list <- list()
for(i in 1:6) {
  dat <- multi_int_fun(i, 365*3, parms, init, amr)[[1]]
  diff_tax_list[[i]] <- data.frame(dat,
                                   "AverageRes" = rowMeans(dat[,4:6]),
                                   "TotInf" = rowSums(dat[,3:6]),
                                   "scen" = paste0("diff", i))
}


test_plot <- melt(diff_tax_list[[6]], id.vars = "time", measure.vars = colnames(diff_tax_list[[6]])[4:6])

ggplot(test_plot, aes(time, value, color = variable)) + geom_line() + theme_bw()

# Obtain the Integrals ----------------------------------------------------

parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5

inf_int = data.frame("Flat_All" = integral(testrun_flat, 3000), 
                     "Dual_12" = integral(dual_list[[1]], 3000), 
                     "Dual_13" = integral(dual_list[[2]], 3000), 
                     "Dual_23" = integral(dual_list[[3]], 3000), 
                     "Single_1" = integral(single_list[[1]], 3000),
                     "Single_2" = integral(single_list[[2]], 3000),
                     "Single_3" = integral(single_list[[3]], 3000),
                     "Diff_1" = integral(diff_tax_list[[1]], 3000),
                     "Diff_2" = integral(diff_tax_list[[2]], 3000),
                     "Diff_3" = integral(diff_tax_list[[3]], 3000),
                     "Diff_4" = integral(diff_tax_list[[4]], 3000),
                     "Diff_5" = integral(diff_tax_list[[5]], 3000),
                     "Diff_6" = integral(diff_tax_list[[6]], 3000))

rownames(inf_int) <- c("Total Infections", "Res_Inf", "WT_Inf", "Average Resistance", "R1", "R2" , "R3")

# Plotting the Matrix -----------------------------------------------------

#All of the plots

OG_melt <- melt(as.matrix(inf_int), measure.vars = colnames(inf_int))
rescale_data <- t(apply(inf_int, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
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

#Limited Number 
lim_OG_melt <- subset(OG_melt, (Var1 == "Total Infections" | Var1 =="Average Resistance"))

ggplot(lim_OG_melt, aes(Var2, Var1)) + theme_bw() +
  geom_tile(aes(fill = rescale)) + 
  facet_grid(Var1 ~ ., scales = "free_y") +
  geom_text(aes(label=value), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "Intervention", expand = c(0, 0))  +   
  scale_y_discrete(name = "Outcome Measure", expand = c(0, 0)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="none")

# Sensitivity Analysis - Effect of Transmission Costs of Resistance -------


# Integrals ---------------------------------------------------------------

diff_tax_list <- list()

for(i in 1:6) {
  diff_tax_list[[i]] <- data.frame(multi_int_fun(i, 365*3, parms, init, amr)[[2]], 
                                   "scen" = paste0("diff", i))
}

#Find The Resistance Lost 
#Baseline

base <- remNA_func(data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms)))
base_totinf <- integral(base, 3000)[1]
base_avgres_int <- integral(base, 3000)[4]

#Lost Resistance Function 
#Supply a dataframe with time and the effective tax over time for class 1, 2 and 3.
#I need to extract the 1-6 data with the 

#Now I need to do this for all 6 scenarios. 
parms_flat <- parms; parms_flat[grep("eff_tax", names(parms), value =T)]  <- 0.5
parms_single1 <- parms; parms_single1[grep("eff_tax1", names(parms), value =T)]  <- 0.5
parms_single2 <- parms; parms_single2[grep("eff_tax2", names(parms), value =T)]  <- 0.5
parms_single3 <- parms; parms_single3[grep("eff_tax3", names(parms), value =T)]  <- 0.5
parms_dual12 <- parms; parms_dual12[grep("eff_tax1|eff_tax2", names(parms), value =T)]  <- 0.5
parms_dual23 <- parms; parms_dual23[grep("eff_tax2|eff_tax3", names(parms), value =T)]  <- 0.5
parms_dual13 <- parms; parms_dual13[grep("eff_tax1|eff_tax2", names(parms), value =T)]  <- 0.5

parm_list <- list(parms_flat,parms_single1,parms_single2,parms_single3,parms_dual12,parms_dual23,parms_dual13)

over_parms <- append(parm_list, diff_tax_list)
list_usage <- lapply(over_parms, usage_fun)

# Relative Reduction Integral Comparison --------------------------------------------

reduc_usage_vec <- sapply(1:length(list_usage), function(x) sum(list_usage[[x]][,6]))
rel_AvgRes <- base_avgres_int - inf_int[4,]
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

# Demonstration of Integrals Relative --------------------------------------

#Plot the Total Infections
#Resistance Integrals
res <- data.frame("time" = base$time,
                  "base_avgres" = rowMeans(base[4:6]),
                  "flat_avgres" = rowMeans(testrun_flat[4:6]))

res$minline = pmin(res[,2], res[,3])
res_df <- melt(res, id.vars=c("time","minline"), variable.name="Scenario", value.name="Resistance", 
               measure.vars = colnames(res)[2:3])

#Total Infection Integrals
totinf <- data.frame("time" = testrun_flat$time,
                     "base_totinf" = rowSums(base[3:6]), 
                     "flat_totinf" = rowSums(testrun_flat[3:6]))

totinf$minline = pmin(totinf[,2], totinf[,3])
totinf_df <- melt(totinf, id.vars=c("time", "minline"), variable.name="Scenario", value.name="TotInf", 
                  measure.vars = colnames(totinf)[2:3])

#Plot Average Resistance
ggplot(res_df, aes(time, Resistance, color = Scenario)) + geom_line(size = 1.2) + theme_bw() + 
  geom_ribbon(aes(ymax=Resistance, ymin=minline), fill = "red", alpha = 0.2) + 
  scale_y_continuous(name = "Average Resistance", expand = c(0, 0), limits = c(0, 0.3))

#Plot Total Infections 
ggplot(totinf_df, aes(time, TotInf, color = Scenario)) + geom_line(size = 1.2) + theme_bw() + 
  geom_ribbon(aes(ymax=TotInf, ymin=minline), fill = "red", alpha = 0.2) + 
  scale_y_continuous(name = "Total Infections", expand = c(0, 0), limits = c(0.9, 1))

#Plot Usage 
usage <- data.frame("time" = testrun_flat$time,
                    "base_usage" = 0.25*3,
                    "usage" = c(rep(0.25*3, 3000),
                                rep((0.25*3)*0.5, 7001)))

usage$minline = pmin(usage[,2], usage[,3])
usage_df <- melt(usage, id.vars=c("time","minline"), variable.name="Scenario", value.name="Usage")

ggplot(usage_df, aes(time, Usage, color = Scenario)) + geom_line(size = 1.2) + theme_bw() + 
  scale_y_continuous(name = "Fraction of Population Treated", limits = c(0,1)) + 
  geom_ribbon(aes(ymax=Usage, ymin=minline), fill = "red", alpha = 0.2)


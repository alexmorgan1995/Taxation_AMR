library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())
setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:10]),
                      sum(data_temp[4:10])),
                      5)
  return(out_vec)
}

# 3 Strain Bonhoeffer ODEs ------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    f1_base <- f_1
    f2_base <- f_2
    f3_base <- f_3
    
    if(t > t_n) {
      f_1 <- f1_base*(1-(eff_tax1_1*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_1*PED2))
      f_3 <- f3_base*(1-(eff_tax3_1*PED3))
    }
    
    if(t > (t_n + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_2*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_2*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_2*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_3*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_3*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_3*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_4*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_4*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_4*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_5*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_5*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_5*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_6*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_6*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_6*PED3)) 
    }

    f_1 <- ifelse(f_1 > 0, f_1, 0)
    f_2 <- ifelse(f_2 > 0, f_2, 0)
    f_3 <- ifelse(f_3 > 0, f_3, 0)
    
    dx = lambda - lambda*x - beta*x*(y_w + y_1 + y_2 + y_3 + y_12 + y_13 + y_23 + y_123) + 
      r_w*y_w + r_1*y_1 + r_2*y_2 + r_3*y_3 + r_12*y_12 + r_23*y_23 + r_13*y_13 + r_123*y_123 +
      h*(1-s)*((f_1+f_2+f_3)*y_w + f_1*(y_2+y_3+y_23)+ f_2*(y_1+y_3+y_13) + f_3*(y_1+y_2+y_12))
    
    dy_w = (beta*x - lambda - r_w - h*(f_1 + f_2 + f_3))*y_w
    
    dy_1 = (beta*x - lambda - r_1 - h*(f_2 + f_3))*y_1 + h*s*f_1*y_w
    dy_2 = (beta*x - lambda - r_2 - h*(f_1 + f_3))*y_2 + h*s*f_2*y_w
    dy_3 = (beta*x - lambda - r_3 - h*(f_2 + f_3))*y_3 + h*s*f_3*y_w
    
    dy_12 = (beta*x - lambda - r_12 - h*f_3)*y_12 + h*s*(f_1*y_2 + f_2*y_1)
    dy_13 = (beta*x - lambda - r_13 - h*f_2)*y_13 + h*s*(f_1*y_3 + f_3*y_1)
    dy_23 = (beta*x - lambda - r_23 - h*f_1)*y_23 + h*s*(f_3*y_2 + f_2*y_3)
    
    dy_123 = (beta*x - lambda - r_123)*y_123 + h*s*(f_1*y_23 + f_2*y_13 + f_3*y_12)
    
    #Calculating the Proportion Integrals
    
    return(list(c(dx, dy_w, dy_1, dy_2, dy_3, dy_12, dy_13, dy_23, dy_123)))
  })
}

# Diff Model --------------------------------------------------------------

multi_int_fun <- function(int_gen, basetax, time_between, parms, init, func){
  
  parms["time_between"] <- time_between
  
  #First Run
  testrun <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))
  values <- tail(testrun, 1)[4:6]
  
  parms[grep("eff_tax1", names(parms), value = T)] <- (basetax*(values[1]/values[2]))
  parms[grep("eff_tax2", names(parms), value = T)] <- (basetax*(values[2]/values[2]))
  parms[grep("eff_tax3", names(parms), value = T)] <- (basetax*(values[3]/values[2]))
  parms[parms < 0] <- 0
  
  testrun <- data.frame(ode(y = init, func = func, times = seq(0, 15000), parms = parms))
  
  #Second Run
  if(int_gen >= 2){
    testrun1 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms))
    values1 <- tail(testrun1, 1)[4:6]
    
    low_char <- names(values1)[which.min(values1)]
    high_char <- names(values1)[which.max(values1)]
    med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]
    
    parms[grep(paste0("eff_tax",substr(low_char, 3, 3)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[low_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char, 3, 3)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[med_char],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char, 3, 3)), names(parms), value = T)[-1]] <- (basetax*(tail(testrun1[high_char],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 15000), parms = parms))
  }
  
  if(int_gen >= 3){
    testrun2 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms))
    values2 <- tail(testrun2, 1)[4:6]
    
    low_char1 <- names(values2)[which.min(values2)]
    high_char1 <- names(values2)[which.max(values2)]
    med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
    
    parms[grep(paste0("eff_tax",substr(low_char1, 3, 3)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[low_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char1, 3, 3)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[med_char1],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char1, 3, 3)), names(parms), value = T)[-c(1,2)]] <- (basetax*(tail(testrun2[high_char1],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 15000), parms = parms))
  }
  
  if(int_gen >= 4){
    testrun3 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms))
    values3 <- tail(testrun3, 1)[4:6]
    
    low_char2 <- names(values3)[which.min(values3)]
    high_char2 <- names(values3)[which.max(values3)]
    med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]
    
    parms[grep(paste0("eff_tax",substr(low_char2, 3, 3)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[low_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char2, 3, 3)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[med_char2],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char2, 3, 3)), names(parms), value = T)[-c(1,2,3)]] <- (basetax*(tail(testrun3[high_char2],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 15000), parms = parms))
  }
  
  if(int_gen >= 5){
    testrun4 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms))
    values4 <- tail(testrun4, 1)[4:6]
    
    low_char3 <- names(values4)[which.min(values4)]
    high_char3 <- names(values4)[which.max(values4)]
    med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
    
    parms[grep(paste0("eff_tax",substr(low_char3, 3, 3)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[low_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(med_char3, 3, 3)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[med_char3],1)/values[2]))
    parms[grep(paste0("eff_tax",substr(high_char3, 3, 3)), names(parms), value = T)[-c(1,2,3,4)]] <- (basetax*(tail(testrun4[high_char3],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 15000), parms = parms))
  }
  
  if(int_gen >= 6){
    testrun5 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms))
    values5 <- tail(testrun5, 1)[4:6]
    
    low_char4 <- names(values5)[which.min(values5)]
    high_char4 <- names(values5)[which.max(values5)]
    med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
    
    parms[grep(paste0("eff_tax", substr(low_char4, 3, 3)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[low_char4],1)/values[2]))
    parms[grep(paste0("eff_tax", substr(med_char4, 3, 3)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[med_char4],1)/values[2]))
    parms[grep(paste0("eff_tax", substr(high_char4, 3, 3)), names(parms), value = T)[-c(1,2,3,4,5)]] <- (basetax*(tail(testrun5[high_char4],1)/values[2]))
    parms[parms < 0] <- 0
    
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 15000), parms = parms))
  }
  return(list(testrun, parms))
}

# Baseline Parms ----------------------------------------------------------

init <- c(x = 0.99, y_w = 1-0.99, y_1 = 0, y_2 = 0, y_3 = 0, y_12 = 0, y_13 = 0, y_23 = 0, y_123 = 0)

parms = c(beta = 2, f_1 = 0.25, f_2 = 0.25, f_3 = 0.25,
          lambda = 0.1, h = 0.01, s = 0.7,
          r_w = 1/12.5, 
          r_1 = 1/12, r_2 = 1/11.85, r_3 = 1/11,
          r_12 = 1/10.5, r_13 = 1/10.4, r_23 = 1/10,
          r_123 = 1/9.75,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 7500, time_between = Inf)

# Run the Models ----------------------------------------------------------

#Flat Tax
parms1 <- parms; parms1[grep("eff_tax", names(parms1), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 15000), parms = parms1))

#Dual Tax
dual_p <- matrix(c(1,1, 2,2, 3,3), ncol = 2, nrow = 3) 
dual_list <- list()
for(i in 1:3) {
  parms1 <- parms
  parms1[grep(paste0("eff_tax", dual_p[i,][1], "|", "eff_tax", dual_p[i,][2]), names(parms1), value = TRUE)] <- 0.5
  dat <- data.frame(ode(y = init, func = amr, times = seq(0, 15000), parms = parms1))
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
  single_list[[i]] <- data.frame(ode(y = init, func = amr, times = seq(0, 15000), parms = parms1),
                                 "scen" = paste0("single_", "eff_tax", i))
}

#Diff Tax
diff_tax_list <- list()

for(i in 1:6) {
  dat <- multi_int_fun(i, 0.5, 365*3, parms, init, amr)[[1]]
  diff_tax_list[[i]] <- data.frame(dat,
                                   "AverageRes" = rowMeans(dat[,4:10]),
                                   "TotInf" = rowSums(dat[,3:10]),
                                   "scen" = paste0("diff", i))
}



# Integrals ---------------------------------------------------------------

inf_int = data.frame("Flat_All" = integral(testrun_flat, 7500), 
                     "Dual_12" = integral(dual_list[[1]], 7500), 
                     "Dual_13" = integral(dual_list[[2]], 7500), 
                     "Dual_23" = integral(dual_list[[3]], 7500), 
                     "Single_1" = integral(single_list[[1]], 7500),
                     "Single_2" = integral(single_list[[2]], 7500),
                     "Single_3" = integral(single_list[[3]], 7500),
                     "Diff_1" = integral(diff_tax_list[[1]], 7500),
                     "Diff_2" = integral(diff_tax_list[[2]], 7500),
                     "Diff_3" = integral(diff_tax_list[[3]], 7500),
                     "Diff_4" = integral(diff_tax_list[[4]], 7500),
                     "Diff_5" = integral(diff_tax_list[[5]], 7500),
                     "Diff_6" = integral(diff_tax_list[[6]], 7500))

rownames(inf_int) <- c("Total Infections", "Average Resistance")

# Plotting the Matrix -----------------------------------------------------

#Limited Number 
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



# Integrals ---------------------------------------------------------------


#Find The Resistance Lost 
#Baseline

base <- data.frame(ode(y = init, func = amr, times = seq(0, 15000), parms = parms))
base_totinf <- integral(base, 7500)[1]
base_avgres_int <- integral(base, 7500)[2]

#Lost Resistance Function 
#Supply a dataframe with time and the effective tax over time for class 1, 2 and 3.
#I need to extract the 1-6 data with the 

#Change in Usage From Baseline - Only for differential taxation models 

usage_fun <- function(parms){
  usage = data.frame("time" = seq(0,7500),
                     "PopUsage1" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax1_", x)]])*0.25, 365*3)),
                                     rep((1-parms[["eff_tax1_6"]])*0.25, 7501-365*3*6)),
                     "PopUsage2" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax2_", x)]])*0.25, 365*3)),
                                     rep((1-parms[["eff_tax2_6"]])*0.25, 7501-365*3*6)),
                     "PopUsage3" = c(sapply(1:6, function(x) rep((1-parms[[paste0("eff_tax3_", x)]])*0.25, 365*3)),
                                     rep((1-parms[["eff_tax3_6"]])*0.25, 7501-365*3*6)))
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

parm_list <- list(parms_flat, parms_single1, parms_single2, parms_single3,
                  parms_dual12, parms_dual23, parms_dual13)

diff_tax_parms <- list()
for(i in 1:6) {
  diff_tax_parms[[i]] <- data.frame(multi_int_fun(i, 0.5, 365*3, parms, init, amr)[[2]], 
                                   "scen" = paste0("diff", i))
}

over_parms <- append(parm_list, diff_tax_parms)
list_usage <- lapply(over_parms, usage_fun)

# Plotting ----------------------------------------------------------------

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
  geom_tile(aes(fill = rescale)) + facet_grid(Var1 ~ ., scales = "free_y") +
  geom_text(aes(label=value), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "Intervention", expand = c(0, 0))  +   
  scale_y_discrete(name = "Outcome Measure", expand = c(0, 0)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="none")

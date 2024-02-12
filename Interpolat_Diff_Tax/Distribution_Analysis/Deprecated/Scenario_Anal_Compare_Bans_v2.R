library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Bans/")


# Calculating Intervention Failure ----------------------------------------

intfail <- function(win_import_change) {
  data <- win_import_change[,c(16:30)]
  prop_vec <- data.frame("intervention" = colnames(data),
                         "prop_inc" = NA)
  for(i in 1:15) {
    prop_1000 <- data[,i]
    if(sum(!is.na(prop_1000)) == 0) {
      prop_vec[i,2] <- NA
    } else{
      prop_vec[i,2] <- length(prop_1000[prop_1000 == -1000])/sum(!is.na(prop_1000)) #out of the interventions which actually run - which ones fail 
    }
  }
  return(prop_vec)
}

# Import in Dataset -------------------------------------------------------

#4 antibiotic classes
win_import_4class <- readRDS("MDR_run_four_tax.RDS")
colnames(win_import_4class)[grep("MR1", colnames(win_import_4class))] <- c("singleMR_inf", "banMR_inf", "singleMR_res" , "banMR_res", 
                                                                           "singleMR_shan","banMR_shan","banMR_avganti")
winintfail_fourclass <- win_import_4class
win_import_4class[win_import_4class == -1000] <- NA
win_import_4class$scen <- "Four Classes"

#Baseline
win_import_base <- readRDS("MDR_run_Base_tax.RDS") 
win_import_base[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                  "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_base$scen <- "Baseline"
win_import_base <- win_import_base[,colnames(win_import_4class)]; winintfail_base <- win_import_base
win_import_base[win_import_base == -1000] <- NA

#Tax 10% of Revenue
win_import_tax10 <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Revisions/MDR_run_Base_tax10.RDS")
win_import_tax10[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                   "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_tax10$scen <- "Tax 10%"
win_import_tax10 <- win_import_tax10[,colnames(win_import_4class)]; winintfail_tax10 <- win_import_tax10
win_import_tax10[win_import_tax10 == -1000] <- NA

#Tax 25% of Revenue
win_import_tax25 <- readRDS("MDR_run_Base_tax25.RDS")
win_import_tax25[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                   "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_tax25$scen <- "Tax 25%"
win_import_tax25 <- win_import_tax25[,colnames(win_import_4class)]; winintfail_tax25 <- win_import_tax25
win_import_tax25[win_import_tax25 == -1000] <- NA

#Tax 75% of Revenue
win_import_tax75 <- readRDS("MDR_run_Base_tax75.RDS")
win_import_tax75[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                   "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_tax75$scen <- "Tax 75%"
win_import_tax75 <- win_import_tax75[,colnames(win_import_4class)]; winintfail_tax75 <- win_import_tax75
win_import_tax75[win_import_tax75 == -1000] <- NA

#Tax 90% of Revenue
win_import_tax90 <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Revisions/MDR_run_Base_tax90.RDS")
win_import_tax90[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                   "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_tax90$scen <- "Tax 90%"
win_import_tax90 <- win_import_tax90[,colnames(win_import_4class)]; winintfail_tax90 <- win_import_tax90
win_import_tax90[win_import_tax90 == -1000] <- NA

#PED of 0.6
win_import_highComp <- readRDS("MDR_run_highComp_tax.RDS")
win_import_highComp[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                      "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_highComp$scen <- "High Cross Elasticity"
win_import_highComp <- win_import_highComp[,colnames(win_import_4class)]; winintfail_highComp <- win_import_highComp
win_import_highComp[win_import_highComp == -1000] <- NA

#PED of 0.4
win_import_lowComp <- readRDS("MDR_run_lowComp_tax.RDS")
win_import_lowComp[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                     "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_lowComp$scen <- "Low Cross Elasticity"
win_import_lowComp <- win_import_lowComp[,colnames(win_import_4class)]; winintfail_lowComp <- win_import_lowComp
win_import_lowComp[win_import_lowComp == -1000] <- NA

#PED Analysis
win_import_PED <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Revisions/MDR_run_PED.RDS")
win_import_PED[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_PED$scen <- "PED Sampling"
win_import_PED <- win_import_PED[,colnames(win_import_4class)]; winintfail_PED <- win_import_PED
win_import_PED[win_import_PED == -1000] <- NA

#Effectiveness threshold of 35%
win_import_35 <- readRDS("MDR_run_35_tax.RDS")
win_import_35[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_35$scen <- "35% Threshold"
win_import_35 <- win_import_35[,colnames(win_import_4class)]; winintfail_tax35 <- win_import_35
win_import_35[win_import_35 == -1000] <- NA


#Effectiveness threshold of 10%
win_import_10 <- readRDS("MDR_run_10_tax.RDS")
win_import_10[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_10$scen <- "10% Threshold"
win_import_10 <- win_import_10[,colnames(win_import_4class)]; winintfail_tax10 <- win_import_10
win_import_10[win_import_10 == -1000] <- NA

#Effectiveness threshold of 5%
win_import_05 <- readRDS("MDR_run_05_tax.RDS")
win_import_05[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_05$scen <- "5% Threshold"
win_import_05 <- win_import_05[,colnames(win_import_4class)]; winintfail_tax05 <- win_import_05
win_import_05[win_import_05 == -1000] <- NA




#Extensive Cattle
win_import_ext_cattle <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Revisions/MDR_cattle_extensive.RDS")
win_import_ext_cattle[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                   "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_ext_cattle$scen <- "Extensive Cattle"
win_import_ext_cattle <- win_import_ext_cattle[,colnames(win_import_4class)]; winintfail_extcattle <- win_import_ext_cattle
win_import_ext_cattle[win_import_ext_cattle == -1000] <- NA

#Extensive Chicken
win_import_ext_chicken <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Revisions/MDR_chick_extense.RDS")
win_import_ext_chicken[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                        "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_ext_chicken$scen <- "Extensive Chicken"
win_import_ext_chicken <- win_import_ext_chicken[,colnames(win_import_4class)]; winintfail_extchick <- win_import_ext_chicken
win_import_ext_chicken[win_import_ext_chicken == -1000] <- NA

#Intensive Cattle
win_import_int_cattle <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Revisions/MDR_cattle_intensive.RDS")
win_import_int_cattle[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                        "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_int_cattle$scen <- "Intensive Cattle"
win_import_int_cattle <- win_import_int_cattle[,colnames(win_import_4class)]; winintfail_intcattle <- win_import_int_cattle
win_import_int_cattle[win_import_int_cattle == -1000] <- NA

#Intensive Chickens
win_import_int_chick <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Revisions/MDR_chick_intensive.RDS")
win_import_int_chick[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                        "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_int_chick$scen <- "Intensive Chicken"
win_import_int_chick <- win_import_int_chick[,colnames(win_import_4class)]; winintfail_intchick <- win_import_int_chick
win_import_int_chick[win_import_int_chick == -1000] <- NA


#2 antibiotic classes
win_import_2class <- readRDS("MDR_run_two_tax.RDS")
win_import_2class[c("singleMR_inf",  "singleMR_res", "singleMR_shan", "singleMR_avganti",
                    "singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                    "banMR_inf",  "banMR_res", "banMR_shan", "banMR_avganti",
                    "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
win_import_2class$scen <- "Two Classes"
win_import_2class <- win_import_2class[,colnames(win_import_4class)]; winintfail_twoclass <- win_import_2class
win_import_2class[win_import_2class == -1000] <- NA


#Intervention Failure
vec_fail <- do.call("rbind", lapply(list(winintfail_base, 
                                         winintfail_tax10, winintfail_tax25, winintfail_tax75, winintfail_tax90,
                                         winintfail_lowComp, winintfail_highComp, winintfail_PED,
                                         winintfail_extcattle, winintfail_intcattle, winintfail_extchick, winintfail_intchick,
                                         winintfail_tax05, winintfail_tax10, winintfail_tax35, 
                                         winintfail_twoclass, winintfail_fourclass), intfail))
#Combine them together
comb_imp <- rbind(win_import_base,
                  win_import_tax10, win_import_tax25, win_import_tax75, win_import_tax90,
                  win_import_lowComp, win_import_highComp, win_import_PED, 
                  winintfail_extcattle, winintfail_extchick, winintfail_intcattle, winintfail_intchick,
                  win_import_05, win_import_10, win_import_05,  win_import_35, 
                  win_import_2class, win_import_4class)

# Altering Data Infections ------------------------------------------------

#Infections
#We want to minimise this - minimise the increase in total infections for each change in usage 
win_inf <- comb_imp[,1:15]

win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

win_inf_trans <- data.frame(win_inf_trans); win_inf_trans$scen <- comb_imp[, 61]

inf_frame <- data.frame(matrix(NA, nrow = 0, ncol = 15))

for(i in unique(win_inf_trans$scen)) {
  data <- win_inf_trans[win_inf_trans$scen == i, 1:15]
  
  prop_win_inf <- data.frame("Infections" = colSums(data)/nrow(data),
                             "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR1)","ST (MR2)",
                                                           "ST (LR)", 
                                                           "DT (1Rd)", "DT (2Rd)",
                                                           "DT (3Rd)", "DT (4Rd)", 
                                                           "DT (5Rd)", "DT (6Rd)",
                                                           "Ban (HR)", "Ban (MR1)","Ban (MR2)", "Ban (LR)")))
  prop_win_inf$Interventions <- factor(prop_win_inf$Interventions, levels = c(prop_win_inf$Interventions))
  prop_win_inf$scen <- i
  inf_frame <- rbind(inf_frame, prop_win_inf)
}

inf_frame$scen <- factor(inf_frame$scen, levels = rev(unique(inf_frame$scen)))

# inf_frame$factor <- c(rep("Baseline", 15), rep("Other", 9*15))
# inf_frame$factor <- factor(inf_frame$factor, levels = unique(inf_frame$factor))

inf_frame$factor <- c(rep("Baseline", 15), rep("Tax", 15*4), rep("PED", 15*3), rep("Produc", 15*4), rep("Thresh", 15*3),  rep("Antibiotics", 15*2))
inf_frame$factor <- factor(inf_frame$factor, levels = unique(inf_frame$factor))

inf_frame$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 17))
inf_frame$Int <- factor(inf_frame$Int, levels = unique(inf_frame$Int))

#Multiplying the Dataframes together
inf_frame$Infections_Scale <- inf_frame$Infections*(1-vec_fail$prop_inc)

inf_scale_vec = c()

for(i in 1:17){
  if(i == 1) {
    reformat <- inf_frame$Infections_Scale[c(1:15)]/sum(inf_frame$Infections_Scale[c(1:15)], na.rm = T)
  } else{
    reformat <- inf_frame$Infections_Scale[c(1:15)+(15*(i-1))]/sum(inf_frame$Infections_Scale[c(1:15)+(15*(i-1))], na.rm = T)
  }
  inf_scale_vec <- append(inf_scale_vec,reformat)
}

inf_frame$Infections_Scale_v1 <- inf_scale_vec

inf_plot <- ggplot(inf_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = Infections_Scale_v1*100)) + scale_fill_distiller(palette ="Blues", direction = 1) +
    facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0)) + ggtitle("Overall Infections (Sensitive + Resistant)") +  
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 0, hjust=0.5)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability of Best Performing Intervention (%)",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Altering Data Resistance ------------------------------------------------

#We want to minimise this - minimise the increase in total infections for each change in usage 
win_res <- comb_imp[,16:30]
#win_inf <- (win_import_pess[,1:10])
#win_inf <- (win_import_opt[,1:10])
win_res[win_res == -1000] <- 1000

win_res_trans <- t(apply(win_res[,1:15], 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

win_res_trans <- data.frame(win_res_trans); win_res_trans$scen <- comb_imp[, 61]

res_frame <- data.frame(matrix(NA, nrow = 0, ncol = 11))

data <- win_res_trans[win_res_trans$scen == "Baseline",1:15]

for(i in unique(win_res_trans$scen)) {
  data <- win_res_trans[win_res_trans$scen == i,1:15]
  prop_win_res <- data.frame("Resistance" = colSums(data, na.rm = T)/nrow(data),
                             "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR1)","ST (MR2)",
                                                           "ST (LR)", 
                                                           "DT (1Rd)", "DT (2Rd)",
                                                           "DT (3Rd)", "DT (4Rd)", 
                                                           "DT (5Rd)", "DT (6Rd)",
                                                           "Ban (HR)", "Ban (MR1)","Ban (MR2)", "Ban (LR)")))
  prop_win_res$Interventions <- factor(prop_win_res$Interventions, levels = c(prop_win_res$Interventions))
  prop_win_res$scen <- i
  res_frame <- rbind(res_frame, prop_win_res)
}

res_frame$scen <- factor(res_frame$scen, levels = rev(unique(res_frame$scen)))
# res_frame$factor <- c(rep("Baseline", 15), rep("Other", 9*15))
# res_frame$factor <- factor(res_frame$factor, levels = unique(res_frame$factor))
res_frame$factor <- c(rep("Baseline", 15), rep("Tax", 15*4), rep("PED", 15*3), rep("Produc", 15*4), rep("Thresh", 15*3),  rep("Antibiotics", 15*2))
res_frame$factor <- factor(res_frame$factor, levels = unique(res_frame$factor))

res_frame$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 17))
res_frame$Int <- factor(res_frame$Int, levels = unique(res_frame$Int))

#Multiplying the Dataframes together
res_frame$Resistance_Scale <- res_frame$Resistance*(1-vec_fail$prop_inc)

res_scale_vec = c()

for(i in 1:17){
  if(i == 1) {
    reformat <- res_frame$Resistance_Scale[c(1:15)]/sum(res_frame$Resistance_Scale[c(1:15)], na.rm = T)
  } else{
    reformat <- res_frame$Resistance_Scale[c(1:15)+(15*(i-1))]/sum(res_frame$Resistance_Scale[c(1:15)+(15*(i-1))], na.rm = T)
  }
  res_scale_vec <- append(res_scale_vec,reformat)
}

res_frame$Resistance_Scale_v1 <- res_scale_vec

res_plot <- ggplot(res_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = Resistance_Scale_v1*100)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0)) +  ggtitle("Average Resistance") +  
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        plot.title = element_text(face = "bold", size = 15),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x =element_text(angle = 0, hjust=0.5)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability of Best Performing Intervention (%)",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Average Antibiotics Available -------------------------------------------

win_avganti <- round((comb_imp[,46:60]), 5)
win_avganti$scen <- comb_imp[, 61]
win_avganti$rowsum <- rowSums(win_avganti[,-16], na.rm = T)
win_avganti <- win_avganti[win_avganti$rowsum!=0,]; scen_remov <- win_avganti[,16]
win_avganti <- win_avganti[,-c(16:17)]

win_avganti_trans <- t(apply(win_avganti, 1, function(x) {
  val = max(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

win_avganti_trans <- data.frame(win_avganti_trans); win_avganti_trans$scen <- scen_remov

avganti_frame <- data.frame(matrix(NA, nrow = 0, ncol = 15))

for(i in unique(win_avganti_trans$scen)) {
  data <- win_avganti_trans[win_avganti_trans$scen == i,1:15]
  
  prop_win_avganti <- data.frame("AverageAnti" = colSums(data)/nrow(data),
                                 "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR1)","ST (MR2)",
                                                               "ST (LR)", 
                                                               "DT (1Rd)", "DT (2Rd)",
                                                               "DT (3Rd)", "DT (4Rd)", 
                                                               "DT (5Rd)", "DT (6Rd)",
                                                               "Ban (HR)", "Ban (MR1)","Ban (MR2)", "Ban (LR)")))
  prop_win_avganti$Interventions <- factor(prop_win_avganti$Interventions, levels = c(prop_win_avganti$Interventions))
  prop_win_avganti$scen <- i
  avganti_frame <- rbind(avganti_frame, prop_win_avganti)
}

avganti_frame$scen <- factor(avganti_frame$scen, levels = rev(unique(avganti_frame$scen)))
# avganti_frame$factor <- c(rep("Baseline", 15), rep("Other", 15*9))
# avganti_frame$factor <- factor(avganti_frame$factor, levels = unique(avganti_frame$factor))
avganti_frame$factor <- c(rep("Baseline", 15), rep("Tax", 15*4), rep("PED", 15*3), rep("Produc", 15*4), rep("Thresh", 15*3),  rep("Antibiotics", 15*2))
avganti_frame$factor <- factor(avganti_frame$factor, levels = unique(avganti_frame$factor))
avganti_frame$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 17))
avganti_frame$Int <- factor(avganti_frame$Int, levels = unique(avganti_frame$Int))

#Multiplying the Dataframes together
avganti_frame$AverageAnti_Scale <- avganti_frame$AverageAnti*(1-vec_fail$prop_inc)

avganti_frame_scale_vec = c()

for(i in 1:17){
  if(i == 1) {
    reformat <- avganti_frame$AverageAnti_Scale[c(1:15)]/sum(avganti_frame$AverageAnti_Scale[c(1:15)], na.rm = T)
  } else{
    reformat <- avganti_frame$AverageAnti_Scale[c(1:15)+(15*(i-1))]/sum(avganti_frame$AverageAnti_Scale[c(1:15)+(15*(i-1))], na.rm = T)
  }
  avganti_frame_scale_vec <- append(avganti_frame_scale_vec,reformat)
}

avganti_frame$avganti_frame_v1 <- avganti_frame_scale_vec

avg_anti_plot <- ggplot(avganti_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = avganti_frame_v1*100)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0))  +  ggtitle("Average Antibiotics Available") +
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom", 
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 0, hjust=0.5)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability of Best Performing Intervention (%)",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Intervention Failure Plot  ----------------------------------------------

vec_fail$scen <- inf_frame$scen

vec_fail$intervention <- rep(c("FT", "ST (HR)", "ST (MR1)","ST (MR2)",
                               "ST (LR)", 
                               "DT (1Rd)", "DT (2Rd)",
                               "DT (3Rd)", "DT (4Rd)", 
                               "DT (5Rd)", "DT (6Rd)",
                               "Ban (HR)", "Ban (MR1)","Ban (MR2)", "Ban (LR)"), 17)

vec_fail$intervention <- factor(vec_fail$intervention, levels = unique(vec_fail$intervention))
# vec_fail$factor <- c(rep("Baseline", 15), rep("Other", 15*9))
# vec_fail$factor <- factor(vec_fail$factor, levels = unique(vec_fail$factor))
vec_fail$factor <- c(rep("Baseline", 15), rep("Tax", 15*4), rep("PED", 15*3), rep("Produc", 15*4), rep("Thresh", 15*3),  rep("Antibiotics", 15*2))
vec_fail$factor <- factor(vec_fail$factor, levels = unique(vec_fail$factor))
vec_fail$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 17))
vec_fail$Int <- factor(vec_fail$Int, levels = unique(vec_fail$Int))

int_failure_p <- ggplot(vec_fail, aes(intervention, scen)) + theme_bw() +
  geom_tile(aes(fill = prop_inc)) + scale_fill_distiller(palette ="Reds", direction = 1)  +
  facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0))  +  ggtitle("Probability of Intervention Failure") +
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom", 
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 0, hjust=0.5)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability of Intervention Failure",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Non-Scaled Plots --------------------------------------------------------

inf_noscale_plot <- ggplot(inf_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = Infections)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0)) + ggtitle("Overall Infections (Sensitive + Resistant)") +  
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 0, hjust=0.5)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

res_plot_noscale_plot <- ggplot(res_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = Resistance)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0)) +  ggtitle("Average Resistance") +  
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        plot.title = element_text(face = "bold", size = 15),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x =element_text(angle = 0, hjust=0.5)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

avg_anti_noscale_plot <- ggplot(avganti_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = AverageAnti)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0))  +  ggtitle("Average Antibiotics Available") +
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom", 
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 0, hjust=0.5)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Combination Plot --------------------------------------------------------

test_scale <- ggarrange(res_plot, inf_plot, avg_anti_plot, ncol = 1, nrow = 3)
test_no_scale <- ggarrange(res_plot_noscale_plot, inf_noscale_plot, avg_anti_noscale_plot, ncol = 1, nrow = 3)

ggsave(test_scale, filename = "scen_compare_pres.png", dpi = 300, width = 13, height = 14, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

ggsave(test_scale, filename = "scen_compare_v1.png", dpi = 300, width = 13, height = 16, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

ggsave(test_no_scale, filename = "scen_compare_NOSCALE_v1.png", dpi = 300, width = 13, height = 16, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

ggsave(int_failure_p, filename = "int_failure_v1.png", dpi = 300, width = 13, height = 6, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

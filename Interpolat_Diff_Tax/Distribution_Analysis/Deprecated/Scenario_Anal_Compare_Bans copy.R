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
      prop_vec[i,2] <- length(prop_1000[prop_1000 == -1000])/sum(!is.na(prop_1000))
    }
  }
  return(prop_vec)
}

# Import in Dataset -------------------------------------------------------

#Baseline
win_import_base <- readRDS("MDR_run_Base_tax.RDS") 
win_import_base[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                  "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_base <- win_import_base
win_import_base[win_import_base == -1000] <- NA
win_import_base$scen <- "Baseline"

#Tax 25% of Revenue
win_import_tax25 <- readRDS("MDR_run_Base_tax25.RDS")
win_import_tax25[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                   "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_tax25 <- win_import_tax25
win_import_tax25[win_import_tax25 == -1000] <- NA
win_import_tax25$scen <- "Tax 25%"

#Tax 75% of Revenue
win_import_tax75 <- readRDS("MDR_run_Base_tax75.RDS")
win_import_tax75[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                   "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_tax75 <- win_import_tax75
win_import_tax75[win_import_tax75 == -1000] <- NA
win_import_tax75$scen <- "Tax 75%"

#PED of 0.6
win_import_highComp <- readRDS("MDR_run_highComp_tax.RDS")
win_import_highComp[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                      "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_highComp <- win_import_highComp
win_import_highComp[win_import_highComp == -1000] <- NA
win_import_highComp$scen <- "High Comp"

#PED of 0.4
win_import_lowComp <- readRDS("MDR_run_lowComp_tax.RDS")
win_import_lowComp[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                     "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_lowComp <- win_import_lowComp
win_import_lowComp[win_import_lowComp == -1000] <- NA
win_import_lowComp$scen <- "Low Comp"

#Effectiveness threshold of 35%
win_import_35 <- readRDS("MDR_run_35_tax.RDS")
win_import_35[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_tax35 <- win_import_35
win_import_35[win_import_35 == -1000] <- NA
win_import_35$scen <- "35% Threshold"

#Effectiveness threshold of 10%
win_import_10 <- readRDS("MDR_run_10_tax.RDS")
win_import_10[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_tax10 <- win_import_10
win_import_10[win_import_10 == -1000] <- NA
win_import_10$scen <- "10% Threshold"

#Effectiveness threshold of 5%
win_import_05 <- readRDS("MDR_run_05_tax.RDS")
win_import_05[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_tax05 <- win_import_05
win_import_05[win_import_05 == -1000] <- NA
win_import_05$scen <- "5% Threshold"

#2 antibiotic classes
win_import_2class <- readRDS("MDR_run_two_tax.RDS")
win_import_2class[c("singleMR_inf",  "singleMR_res", "singleMR_shan", "singleMR_avganti",
                    "singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti",
                    "banMR_inf",  "banMR_res", "banMR_shan", "banMR_avganti",
                    "banMR2_inf",  "banMR2_res", "banMR2_shan", "banMR2_avganti")] <- NA
winintfail_twoclass <- win_import_2class
win_import_2class[win_import_2class == -1000] <- NA
win_import_2class$scen <- "Two Classes"

#4 antibiotic classes
win_import_4class <- readRDS("MDR_run_four_tax.RDS")
colnames(win_import_4class)[grep("MR1", colnames(win_import_4class))] <- c("singleMR_inf", "banMR_inf", "singleMR_res" , "banMR_res", 
                                                                           "singleMR_shan","banMR_shan","banMR_avganti")
winintfail_fourclass <- win_import_4class
win_import_4class[win_import_4class == -1000] <- NA
win_import_4class$scen <- "Four Classes"


#The MR1 in the avganti for single MR intervention is already without the MR1

win_import_base <- win_import_base[,colnames(win_import_4class)]; fail_base <- intfail(winintfail_base)
win_import_tax25 <- win_import_tax25[,colnames(win_import_4class)]; fail_tax25 <- intfail(winintfail_tax25)
win_import_tax75 <- win_import_tax75[,colnames(win_import_4class)]; fail_tax75 <- intfail(winintfail_tax75)
win_import_highComp <- win_import_highComp[,colnames(win_import_4class)]; fail_highComp <- intfail(winintfail_highComp)
win_import_lowComp <- win_import_lowComp[,colnames(win_import_4class)]; fail_lowComp <- intfail(win_import_lowComp)
win_import_35 <- win_import_35[,colnames(win_import_4class)]; fail_35 <- intfail(winintfail_tax35)
win_import_10 <- win_import_10[,colnames(win_import_4class)]; fail_10 <- intfail(winintfail_tax10)
win_import_05 <- win_import_05[,colnames(win_import_4class)]; fail_05 <- intfail(winintfail_tax05)
win_import_2class <- win_import_2class[,colnames(win_import_4class)]; fail_2class <- intfail(winintfail_twoclass)
fail_4class <- intfail(winintfail_fourclass)


intfail(win_import_lowComp)

data <- win_import_lowComp[,c(16:30)]
prop_vec <- data.frame("intervention" = colnames(data),
                       "prop_inc" = NA)
for(i in 1:15) {
  prop_1000 <- data[,i]
  if(sum(!is.na(prop_1000)) == 0) {
    prop_vec[i,2] <- NA
  } else{
    prop_vec[i,2] <- length(prop_1000[prop_1000 == -1000])/sum(!is.na(prop_1000))
  }
}


vec_fail <- do.call("rbind", lapply(list(win_import_base, win_import_tax25, win_import_tax75, win_import_lowComp, win_import_highComp,
                                         win_import_35, win_import_10, win_import_05, win_import_2class, win_import_4class), intfail))
#Combine them together
comb_imp <- rbind(win_import_base,
                  win_import_tax25, win_import_tax75,
                  win_import_lowComp, win_import_highComp, 
                  win_import_35, win_import_10, win_import_05, 
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
inf_frame$factor <- c(rep("Baseline", 15), rep("Other", 9*15))
inf_frame$factor <- factor(inf_frame$factor, levels = unique(inf_frame$factor))
inf_frame$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 10))
inf_frame$Int <- factor(inf_frame$Int, levels = unique(inf_frame$Int))

#Multiplying the Dataframes together
inf_frame$Infections_Scale <- inf_frame$Infections*vec_fail

for(i in 1:10){
  reformat <- inf_frame$Infections_Scale[c(1:15)+(15*i-1)]
  print(reformat)
}



inf_plot <- ggplot(inf_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = Infections)) + scale_fill_distiller(palette ="Blues", direction = 1) +
    facet_grid(factor~ Int, scales  = "free", space = "free") +
  scale_x_discrete(name = "", expand = c(0, 0)) + ggtitle("Total Infections") +  
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

for(i in unique(win_res_trans$scen)) {
  data <- win_res_trans[win_res_trans$scen == i,1:15]
  
  prop_win_res <- data.frame("Resistance" = colSums(data)/nrow(data),
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
res_frame$factor <- c(rep("Baseline", 15), rep("Other", 6*15))
res_frame$factor <- factor(res_frame$factor, levels = unique(res_frame$factor))
res_frame$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 7))
res_frame$Int <- factor(res_frame$Int, levels = unique(res_frame$Int))

res_plot <- ggplot(res_frame, aes(Interventions, scen)) + theme_bw() +
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
avganti_frame$factor <- c(rep("Baseline", 15), rep("Other", 15*6))
avganti_frame$factor <- factor(avganti_frame$factor, levels = unique(avganti_frame$factor))
avganti_frame$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 7))
avganti_frame$Int <- factor(avganti_frame$Int, levels = unique(avganti_frame$Int))

avg_anti_plot <- ggplot(avganti_frame, aes(Interventions, scen)) + theme_bw() +
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

# Intervention Failure Plot  ----------------------------------------------

int_failure <- comb_imp[,c(16:30,61)]

prop_vec <- data.frame("intervention" = as.factor(rep(colnames(int_failure)[-16], 7)),
                       "prop_inc" = NA,
                       "scen" = NA)

for(z in 1:7) {
  win_import <- int_failure[int_failure$scen == unique(int_failure$scen)[z],-16]
  for(i in 1:15) {
    prop_1000 <- win_import[,i]
    prop_vec[i+(15*(z-1)),2] <- length(prop_1000[prop_1000 == -1000])/sum(!is.na(prop_1000))
    prop_vec[i+(15*(z-1)),3] <- c("Baseline", "Case Study", "Severe", "10% Thresh", "5% Thresh", "Two Classes", "Four Classes")[z]
  }
}

prop_vec$intervention <- rep(c("FT", "ST (HR)", "ST (MR1)","ST (MR2)",
                               "ST (LR)", 
                               "DT (1Rd)", "DT (2Rd)",
                               "DT (3Rd)", "DT (4Rd)", 
                               "DT (5Rd)", "DT (6Rd)",
                               "Ban (HR)", "Ban (MR1)","Ban (MR2)", "Ban (LR)"), 7)

prop_vec$intervention <- factor(prop_vec$intervention, levels = unique(prop_vec$intervention))
prop_vec$factor <- c(rep("Baseline", 15), rep("Other", 15*6))
prop_vec$factor <- factor(prop_vec$factor, levels = unique(prop_vec$factor))
prop_vec$Int <- as.factor(rep(c(rep("Tax", 11),rep("Bans", 4)), 7))
prop_vec$Int <- factor(prop_vec$Int, levels = unique(prop_vec$Int))

int_failure_p <- ggplot(prop_vec, aes(intervention, scen)) + theme_bw() +
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

# Combination Plot --------------------------------------------------------

test <- ggarrange(res_plot, inf_plot,avg_anti_plot, ncol = 1, nrow = 3)

ggsave(test, filename = "scen_compare.png", dpi = 300, width = 13, height = 12, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

ggsave(int_failure_p, filename = "int_failure.png", dpi = 300, width = 13, height = 6, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

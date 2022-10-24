library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output")

# Import in Dataset -------------------------------------------------------

win_import_base <- readRDS("MDR_run_interpol.RDS")
win_import_base[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti")] <- NA
win_import_base$scen <- "Baseline"

win_import_25 <- readRDS("MDR_run_interpol_25.RDS")
win_import_25[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti")] <- NA
win_import_25$scen <- "25% Threshold"

win_import_75 <- readRDS("MDR_run_interpol_75.RDS")
win_import_75[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti")] <- NA
win_import_75$scen <- "75% Threshold"

win_import_2class <- readRDS("MDR_run_two.RDS")
win_import_2class[c("singleMR_inf",  "singleMR_res", "singleMR_shan", "singleMR_avganti",
                    "singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti")] <- NA
win_import_2class$scen <- "Two Classes"

win_import_4class <- readRDS("MDR_run_four.RDS")
win_import_4class$scen <- "Four Classes"
colnames(win_import_4class)[grep("MR1", colnames(win_import_4class))] <- c("singleMR_inf",  "singleMR_res" , "singleMR_shan")

win_import_realPED <- readRDS("MDR_run_interpol_realPED.RDS")
win_import_realPED[c("singleMR2_inf",  "singleMR2_res", "singleMR2_shan", "singleMR2_avganti")] <- NA
win_import_realPED$scen <- "Real PED"


win_import_base <- win_import_base[,colnames(win_import_4class)]
win_import_25 <- win_import_25[,colnames(win_import_4class)]
win_import_75 <- win_import_75[,colnames(win_import_4class)]
win_import_2class <- win_import_2class[,colnames(win_import_4class)]
win_import_realPED <- win_import_realPED[,colnames(win_import_4class)]

comb_imp <- rbind(win_import_base, win_import_25, win_import_75, win_import_realPED, win_import_2class, win_import_4class)

# Altering Data Infections ------------------------------------------------

#Infections
#We want to minimise this - minimise the increase in total infections for each change in usage 
win_inf <- comb_imp[,1:11]
#win_inf <- (win_import_pess[,1:10])
#win_inf <- (win_import_opt[,1:10])

win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))


win_inf_trans <- data.frame(win_inf_trans); win_inf_trans$scen <- comb_imp[, 45]

inf_frame <- data.frame(matrix(NA, nrow = 0, ncol = 11))

for(i in unique(win_inf_trans$scen)) {
  data <- win_inf_trans[win_inf_trans$scen == i,1:11]
  
  prop_win_inf <- data.frame("Infections" = colSums(data)/nrow(data),
                             "Interventions" = as.factor(c("Flat", "Single (HR)", "Single (MR)","Single (MR2)",
                                                           "Single (LR)", 
                                                           "Diff (1 Rd)", "Diff (2 Rd)",
                                                           "Diff (3 Rd)", "Diff (4 Rd)", 
                                                           "Diff (5 Rd)", "Diff (6 Rd)")))
  prop_win_inf$Interventions <- factor(prop_win_inf$Interventions, levels = c(prop_win_inf$Interventions))
  prop_win_inf$scen <- i
  inf_frame <- rbind(inf_frame, prop_win_inf)
}

inf_frame$scen <- factor(inf_frame$scen, levels = rev(unique(inf_frame$scen)))

inf_plot <- ggplot(inf_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = Infections)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0)) + ggtitle("Total Infections") +  
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Altering Data Resistance ------------------------------------------------

#Infections
#We want to minimise this - minimise the increase in total infections for each change in usage 
win_res <- comb_imp[,12:22]
#win_inf <- (win_import_pess[,1:10])
#win_inf <- (win_import_opt[,1:10])

win_res_trans <- t(apply(win_res[,1:11], 1, function(x) {
  val = max(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

win_res_trans <- data.frame(win_res_trans); win_res_trans$scen <- comb_imp[, 45]

res_frame <- data.frame(matrix(NA, nrow = 0, ncol = 11))

for(i in unique(win_res_trans$scen)) {
  data <- win_res_trans[win_res_trans$scen == i,1:11]
  
  prop_win_res <- data.frame("Resistance" = colSums(data)/nrow(data),
                             "Interventions" = as.factor(c("Flat", "Single (HR)", "Single (MR)","Single (MR2)",
                                                           "Single (LR)", 
                                                           "Diff (1 Rd)", "Diff (2 Rd)",
                                                           "Diff (3 Rd)", "Diff (4 Rd)", 
                                                           "Diff (5 Rd)", "Diff (6 Rd)")))
  prop_win_res$Interventions <- factor(prop_win_res$Interventions, levels = c(prop_win_res$Interventions))
  prop_win_res$scen <- i
  res_frame <- rbind(res_frame, prop_win_res)
}

res_frame$scen <- factor(res_frame$scen, levels = rev(unique(res_frame$scen)))

res_plot <- ggplot(res_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = Resistance)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0)) +  ggtitle("Average Resistance") +  
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        plot.title = element_text(face = "bold", size = 15),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x =element_text(angle = 45, hjust=1)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Average Antibiotics Available -------------------------------------------

win_avganti <- round((comb_imp[,34:44]), 5)
win_avganti$scen <- comb_imp[, 45]
win_avganti$rowsum <- rowSums(win_avganti[,-12], na.rm = T)
win_avganti <- win_avganti[win_avganti$rowsum!=0,]; scen_remov <- win_avganti[,12]
win_avganti <- win_avganti[,-c(12:13)]

win_avganti_trans <- t(apply(win_avganti, 1, function(x) {
  val = max(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

win_avganti_trans <- data.frame(win_avganti_trans); win_avganti_trans$scen <- scen_remov

avganti_frame <- data.frame(matrix(NA, nrow = 0, ncol = 11))

for(i in unique(win_avganti_trans$scen)) {
  data <- win_avganti_trans[win_avganti_trans$scen == i,1:11]
  
  prop_win_avganti <- data.frame("AverageAnti" = colSums(data)/nrow(data),
                                 "Interventions" = as.factor(c("Flat", "Single (HR)", "Single (MR)","Single (MR2)",
                                                               "Single (LR)", 
                                                               "Diff (1 Rd)", "Diff (2 Rd)",
                                                               "Diff (3 Rd)", "Diff (4 Rd)", 
                                                               "Diff (5 Rd)", "Diff (6 Rd)")))
  prop_win_avganti$Interventions <- factor(prop_win_avganti$Interventions, levels = c(prop_win_avganti$Interventions))
  prop_win_avganti$scen <- i
  avganti_frame <- rbind(avganti_frame, prop_win_avganti)
}

avganti_frame$scen <- factor(avganti_frame$scen, levels = rev(unique(avganti_frame$scen)))

avg_anti_plot <- ggplot(avganti_frame, aes(Interventions, scen)) + theme_bw() +
  geom_tile(aes(fill = AverageAnti)) + scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +  ggtitle("Average Antibiotics Available") +
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom", 
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 10,
                               barheight = 1)) 

# Combination Plot --------------------------------------------------------


ggarrange(res_plot, inf_plot,avg_anti_plot, ncol = 1, nrow = 3)

c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)","Single Tax (MR2)",
  "Single Tax (LR)", 
  "Diff Tax (1 Round)", "Diff Tax (2 Round)",
  "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
  "Diff Tax (5 Round)", "Diff Tax (6 Round)")

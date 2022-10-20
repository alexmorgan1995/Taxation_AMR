library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
library("confintr")
rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Formalised_Analysis/Model_Output/Comparison_Sens")

# Import in Dataset -------------------------------------------------------

win_import_base <- readRDS("MDR_run_v5.RDS")
win_import_075 <- readRDS("MDR_run_v5_075.RDS")
win_import_025 <- readRDS("MDR_run_v5_025.RDS")

# Altering Data Shannon's Index -------------------------------------------

#
baseline_shan <- data.frame(round((win_import_base[,21:30]), 5))
base.data.frame_shan <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:10) {
  base.data.frame_shan <- rbind(base.data.frame_shan, 
                                data.frame("variable" = colnames(baseline_shan)[i],
                                            "mean" = mean(as.numeric(baseline_shan[,i]), na.rm = T),
                                           "lowerCI" = ci_mean(as.numeric(baseline_shan[,i]))[[2]][1],
                                           "upperCI" = ci_mean(as.numeric(baseline_shan[,i]))[[2]][2],
                                           "group" = "baseline"))
}

#
win_075_shan <- data.frame(round((win_import_075[,21:30]), 5))
win_075.data.frame_shan <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:10) {
  win_075.data.frame_shan <- rbind(win_075.data.frame_shan, 
                                   data.frame("variable" = colnames(baseline_shan)[i],
                                              "mean" = mean(as.numeric(win_075_shan[,i]), na.rm = T),
                                              "lowerCI" = ci_mean(as.numeric(win_075_shan[,i]))[[2]][1],
                                              "upperCI" = ci_mean(as.numeric(win_075_shan[,i]))[[2]][2],
                                              "group" = "win_075"))
}

#
win_025_shan <- data.frame(round((win_import_025[,21:30]), 5))
win_025.data.frame_shan <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:10) {
  win_025.data.frame_shan <- rbind(win_025.data.frame_shan, 
                                   data.frame("variable" = colnames(baseline_shan)[i],
                                              "mean" = mean(as.numeric(win_025_shan[,i]), na.rm = T),
                                              "lowerCI" = ci_mean(as.numeric(win_025_shan[,i]))[[2]][1],
                                              "upperCI" = ci_mean(as.numeric(win_025_shan[,i]))[[2]][2],
                                              "group" = "win_025"))
}

comb_shan <- rbind(base.data.frame_shan, win_075.data.frame_shan, win_025.data.frame_shan)

# Average Antibiotics Available -------------------------------------------

#
baseline_avganti <- data.frame(round((win_import_base[,31:40]), 5))
base.data.frame_avganti <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:10) {
  base.data.frame_avganti <- rbind(base.data.frame_avganti, 
                                data.frame("variable" = colnames(baseline_avganti)[i],
                                           "mean" = mean(as.numeric(baseline_avganti[,i]), na.rm = T),
                                           "lowerCI" = ci_mean(as.numeric(baseline_avganti[,i]))[[2]][1],
                                           "upperCI" = ci_mean(as.numeric(baseline_avganti[,i]))[[2]][2],
                                           "group" = "baseline"))
}

#
win_075_avganti <- data.frame(round((win_import_075[,31:40]), 5))
win_075.data.frame_avganti <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:10) {
  win_075.data.frame_avganti <- rbind(win_075.data.frame_avganti, 
                                   data.frame("variable" = colnames(baseline_avganti)[i],
                                              "mean" = mean(as.numeric(win_075_avganti[,i]), na.rm = T),
                                              "lowerCI" = ci_mean(as.numeric(win_075_avganti[,i]))[[2]][1],
                                              "upperCI" = ci_mean(as.numeric(win_075_avganti[,i]))[[2]][2],
                                              "group" = "win_075"))
}

#
win_025_avganti <- data.frame(round((win_import_025[,21:30]), 5))
win_025.data.frame_avganti <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:10) {
  win_025.data.frame_avganti <- rbind(win_025.data.frame_avganti, 
                                   data.frame("variable" = colnames(baseline_avganti)[i],
                                              "mean" = mean(as.numeric(win_025_avganti[,i]), na.rm = T),
                                              "lowerCI" = ci_mean(as.numeric(win_025_avganti[,i]))[[2]][1],
                                              "upperCI" = ci_mean(as.numeric(win_025_avganti[,i]))[[2]][2],
                                              "group" = "win_025"))
}

comb_avganti <- rbind(base.data.frame_avganti, win_075.data.frame_avganti, win_025.data.frame_avganti)


# Plotting the Data ------------------------------------------------------

ggplot(comb_shan, aes(variable, mean)) + theme_bw() +
  geom_point(aes(color = group), position = position_dodge(0.3))


ggplot(comb_avganti, aes(variable, mean)) +
  geom_point(aes(color = group), position = position_dodge(0.3))

library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Formalised_Analysis/Model_Output/Comparison_Sens")

# Import in Dataset -------------------------------------------------------

win_import_base <- readRDS("MDR_run_four.RDS")
win_import_075 <- readRDS("MDR_run_v5_075.RDS")
win_import_025 <- readRDS("MDR_run_v5_025.RDS")

# Altering Data Shannon's Index -------------------------------------------

#
baseline_shan <- data.frame(round((win_import_base[,21:30]), 5))
base.data.frame_shan <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:11) {
  base.data.frame_shan <- rbind(base.data.frame_shan, 
                                data.frame("variable" = colnames(baseline_shan)[i],
                                            "value" = as.numeric(baseline_shan[i,]),
                                           "group" = "baseline"))
}

#
win_075_shan <- data.frame(round((win_import_075[,21:30]), 5))
win_075.data.frame_shan <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:11) {
  win_075.data.frame_shan <- rbind(win_075.data.frame_shan, 
                                   data.frame("variable" = colnames(win_075_shan)[i],
                                              "value" = as.numeric(win_075_shan[i,]),
                                              "group" = "win_075"))
}

#
win_025_shan <- data.frame(round((win_import_025[,21:30]), 5))
win_025.data.frame_shan <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:11) {
  win_025.data.frame_shan <- rbind(win_025.data.frame_shan, 
                                   data.frame("variable" = colnames(win_025_shan)[i],
                                              "value" = as.numeric(win_025_shan[i,]),
                                              "group" = "win_025"))
}


comb_shan <- rbind(base.data.frame_shan, win_075.data.frame_shan, win_025.data.frame_shan)

# Altering Data Shannon's Index -------------------------------------------

#
baseline <- data.frame(round((win_import_base[,23:33]), 5))
base.data.frame <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:11) {
  base.data.frame <- rbind(base.data.frame, 
                           data.frame("variable" = colnames(baseline)[i],
                                      "value" = as.numeric(baseline[i,]),
                                      "group" = "baseline"))
}

#
win_075 <- data.frame(round((win_import_075[,23:33]), 5))
win_075.data.frame <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:11) {
  win_075.data.frame <- rbind(win_075.data.frame, 
                              data.frame("variable" = colnames(win_075)[i],
                                         "value" = as.numeric(win_075[i,]),
                                         "group" = "win_075"))
}

#
win_025 <- data.frame(round((win_import_025[,23:33]), 5))
win_025.data.frame <- data.frame(matrix(ncol = 2, nrow = 0))

for(i in 1:11) {
  win_025.data.frame <- rbind(win_025.data.frame, 
                              data.frame("variable" = colnames(win_025)[i],
                                         "value" = as.numeric(win_025[i,]),
                                         "group" = "win_025"))
}


# Combine the Dataframe ---------------------------------------------------






win_075 <- data.frame("Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                                                     "Single Tax (LR)", 
                                                     "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                     "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                     "Diff Tax (5 Round)", "Diff Tax (6 Round)")),
                       round((win_import_075[,23:33]), 5))

win_075 <- data.frame(melt(win_075, id.vars = "Interventions", measure.vars = colnames(win_075)[-1]),
                       group = "win_075")


win_025 <- data.frame("Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                                                     "Single Tax (LR)", 
                                                     "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                     "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                     "Diff Tax (5 Round)", "Diff Tax (6 Round)")),
                       round((win_import_025[,23:33]), 5))

win_025 <- data.frame(melt(win_025, id.vars = "Interventions", measure.vars = colnames(win_025)[-1]),
                       group = "win_025")





win_shan <- data.frame("Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                                                     "Single Tax (LR)", 
                                                     "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                     "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                     "Diff Tax (5 Round)", "Diff Tax (6 Round)")),
                       "baseline" = round((win_import_base[,23:33]), 5),
                       "thresh_075" = round((win_import_075[,23:33]), 5),
                       "thresh_025" = round((win_import_025[,23:33]), 5))

melt_shan <- melt(win_shan, id.vars = "Interventions", measure.vars = colnames(win_shan)[-1])

#prop_win_shan$Interventions <- factor(prop_win_shan$Interventions, levels = c(prop_win_shan$Interventions))

# Average Antibiotics Available -------------------------------------------

win_avganti <- data.frame("Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                                                     "Single Tax (LR)", 
                                                     "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                     "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                     "Diff Tax (5 Round)", "Diff Tax (6 Round)")),
                       "baseline" = round((win_import[,34:44]), 5),
                       "thresh_075" = round((win_import_075[,34:44]), 5),
                       "thresh_025" = round((win_import_025[,34:44]), 5))

melt_avganti <- melt(win_avganti, id.vars = "Interventions", measure.vars = colnames(win_avganti)[-1])

#prop_win_avganti$Interventions <- factor(prop_win_avganti$Interventions, levels = c(prop_win_avganti$Interventions))

# Combining All Together --------------------------------------------------


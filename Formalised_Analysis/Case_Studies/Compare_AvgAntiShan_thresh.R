library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Formalised_Analysis/Model_Output/Comparison_Sens")

# Import in Dataset -------------------------------------------------------

win_import_base <- readRDS("MDR_run_four.RDS")
win_import_075 <- readRDS("MDR_run_v5_075.RDS")
win_import_025 <- readRDS("MDR_run_v5_025.RDS")

for(i in seq_along(win_import)) {
  win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
}

# Altering Data Shannon's Index -------------------------------------------

win_shan <- data.frame("Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                                                     "Single Tax (LR)", 
                                                     "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                     "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                     "Diff Tax (5 Round)", "Diff Tax (6 Round)")),
                       "baseline" = round((win_import[,23:33]), 5),
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


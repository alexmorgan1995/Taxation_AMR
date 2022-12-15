library("deSolve"); library("parallel"); library("lhs"); library("dplyr")
rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Bans")

# Import in the Data ------------------------------------------------------

win_import_base <- readRDS("MDR_run_ban.RDS")
win_import_parms <- readRDS("MDR_run_parms_ban.RDS")

# Isolate the Average Antibiotics Available -------------------------------

avg_anti <- win_import_base[,40:52]
avg_anti_trim <- avg_anti[rowSums(avg_anti[,10] > 2) > 0, ]

cleanColls <- function(x) {
  x <- x[-(26:28)]
}
avg_anti[avg_anti[,10] > 2,10]
df <- do.call(rbind.data.frame, lapply(win_import_parms, cleanColls))[avg_anti[,10] > 2,-(24:26)]
df$data <- avg_anti[avg_anti[,10] > 2,10]

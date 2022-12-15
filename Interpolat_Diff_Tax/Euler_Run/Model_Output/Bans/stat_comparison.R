library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Bans")

# Import in Dataset -------------------------------------------------------

win_import_change <- readRDS("MDR_run_ban.RDS"); win_import <- win_import_change
win_import_change <- readRDS("MDR_run_ban_75.RDS"); win_import <- win_import_change
win_import_change <- readRDS("MDR_run_ban_25.RDS"); win_import <- win_import_change
win_import_change <- readRDS("MDR_run_ban_bias_PED.RDS"); win_import <- win_import_change
win_import_change <- readRDS("MDR_run_ban_realPED.RDS"); win_import <- win_import_change


for(i in seq_along(win_import)) {
  win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
}

win_import[win_import == -1000] <- NA

win_inf <- (win_import[,1:13])
win_res <- (win_import[,14:26])

#Res

win_res_trans <- t(apply(win_res, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = colSums(win_res_trans, na.rm = T)/nrow(win_res_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (MR)", "Ban (LR)")))

#Infection

win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("Infection" = colSums(win_inf_trans, na.rm = T)/nrow(win_inf_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (MR)",  "Ban (LR)")))

#Plot

m_res <- melt(win_res, measure.vars = colnames(win_res))
m_inf <- melt(win_inf, measure.vars = colnames(win_inf))

test_stat_res <- pairwise.wilcox.test(m_res$value, m_res$variable,
                                      p.adjust.method = "bonferroni")

test_stat_inf <- pairwise.wilcox.test(m_inf$value, m_inf$variable,
                                      p.adjust.method = "bonferroni")

# Two Antibiotics ---------------------------------------------------------

win_import_change <- readRDS("MDR_run_two.RDS"); win_import <- win_import_change

for(i in seq_along(win_import)) {
  win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
}

win_import[win_import == -1000] <- NA

win_inf <- (win_import[,1:11])
win_res <- (win_import[,12:22])

#Res
win_res_trans <- t(apply(win_res, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = colSums(win_res_trans, na.rm = T)/nrow(win_res_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", 
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (LR)")))

#Infection
win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("Infection" = colSums(win_inf_trans, na.rm = T)/nrow(win_inf_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", 
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (LR)")))

#Plot
m_res <- melt(win_res, measure.vars = colnames(win_res))
m_inf <- melt(win_inf, measure.vars = colnames(win_inf))

test_stat_res <- pairwise.wilcox.test(m_res$value, m_res$variable,
                                      p.adjust.method = "bonferroni")

test_stat_inf <- pairwise.wilcox.test(m_inf$value, m_inf$variable,
                                      p.adjust.method = "bonferroni")

# Four Antibiotics --------------------------------------------------------

win_import_change <- readRDS("MDR_run_four.RDS"); win_import <- win_import_change

for(i in seq_along(win_import)) {
  win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
}

win_import[win_import == -1000] <- NA

win_inf <- (win_import[,1:15])
win_res <- (win_import[,16:30])

#Res

win_res_trans <- t(apply(win_res, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = colSums(win_res_trans, na.rm = T)/nrow(win_res_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR1)","ST (MR2)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (MR1)", "Ban (MR2)", "Ban (LR)")))

#Infection

win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("Infection" = colSums(win_inf_trans, na.rm = T)/nrow(win_inf_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR1)","ST (MR2)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (MR1)", "Ban (MR2)", "Ban (LR)")))

#Plot

m_res <- melt(win_res, measure.vars = colnames(win_res))
m_inf <- melt(win_inf, measure.vars = colnames(win_inf))

test_stat_res <- pairwise.wilcox.test(m_res$value, m_res$variable,
                                      p.adjust.method = "bonferroni")

test_stat_inf <- pairwise.wilcox.test(m_inf$value, m_inf$variable,
                                      p.adjust.method = "bonferroni")
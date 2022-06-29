library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/Theoretical_Analysis/Formalised_Anal/Sens_Anal/")

# Import in Dataset -------------------------------------------------------

win_import <- readRDS("comb_data.RDS")
for(i in seq_along(win_import)) win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
# Altering Data -----------------------------------------------------------

#Infections
#We want to minimise this - minimise the increase in total infections for each change in usage 
win_inf <- (win_import[,1:10])

win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("prop_win" = colSums(win_inf_trans)/nrow(win_inf_trans),
                           "interventions" = as.factor(colnames(win_inf_trans)))
prop_win_inf$interventions <- factor(prop_win_inf$interventions, levels = c(prop_win_inf$interventions ))

#Resistance

win_res <- (win_import[,11:20])

win_res_trans <- t(apply(win_res, 1, function(x) {
  val = max(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("prop_win" = colSums(win_res_trans)/nrow(win_res_trans),
                           "interventions" = as.factor(colnames(win_res_trans)))
prop_win_res$interventions <- factor(prop_win_res$interventions, levels = c(prop_win_res$interventions ))

# Plotting win Probabilities ----------------------------------------------

p_inf <- ggplot(prop_win_inf, aes(y = prop_win, x = as.factor(interventions))) + geom_bar(stat="identity")  + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(y ="Probability of Minimising Infection Increase", x = "Interventions", fill = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))

p_res <- ggplot(prop_win_res, aes(y = prop_win, x = as.factor(interventions))) + geom_bar(stat="identity") + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(y ="Probability of Maximising Resistance Decrease", x = "Interventions", fill = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))

ggarrange(p_inf, p_res, ncol = 1, nrow = 2)

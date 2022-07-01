library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/Theoretical_Analysis/Formalised_Anal/Sens_Anal/PED_Scen/")

# Import in Dataset -------------------------------------------------------

win_import <- readRDS("comb_data.RDS")
win_import_pess <- readRDS("comb_data_pess.RDS")
win_import_opt <- readRDS("comb_data_opt.RDS")

for(i in seq_along(win_import)) win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
for(i in seq_along(win_import_pess)) win_import_pess[[i]] <- as(win_import_pess[[i]], class(win_import_pess[[i]][[1]]))
for(i in seq_along(win_import_opt)) win_import_opt[[i]] <- as(win_import_opt[[i]], class(win_import_opt[[i]][[1]]))

# Altering Data -----------------------------------------------------------

#Infections
#We want to minimise this - minimise the increase in total infections for each change in usage 
win_inf <- (win_import[,1:10])
win_inf <- (win_import_pess[,1:10])
win_inf <- (win_import_opt[,1:10])



win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("Infections" = colSums(win_inf_trans)/nrow(win_inf_trans),
                           "Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                                                         "Single Tax (LR)", 
                                                         "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                         "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                         "Diff Tax (5 Round)", "Diff Tax (6 Round)")))
prop_win_inf$Interventions <- factor(prop_win_inf$Interventions, levels = c(prop_win_inf$Interventions ))

#Resistance

win_res <- (win_import[,11:20])
win_res <- (win_import_pess[,11:20])
win_res <- (win_import_opt[,11:20])



win_res_trans <- t(apply(win_res, 1, function(x) {
  val = max(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = colSums(win_res_trans)/nrow(win_res_trans),
                           "Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                                                         "Single Tax (LR)", 
                                                         "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                         "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                         "Diff Tax (5 Round)", "Diff Tax (6 Round)")))
prop_win_res$Interventions <- factor(prop_win_res$Interventions, levels = c(prop_win_res$Interventions ))

#Combine the Two Together 

combdata <- prop_win_res; combdata$Infections <- prop_win_inf$Infections

melt_combdata <- melt(combdata, id.vars = "Interventions", measure.vars = c("Resistance", "Infections"))
melt_combdata$Interventions <- factor(melt_combdata$Interventions, levels = c(prop_win_res$Interventions ))

# Plotting win Probabilities ----------------------------------------------

p_inf <- ggplot(prop_win_inf, aes(y = Infections, x = as.factor(Interventions))) + geom_bar(stat="identity")  + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(y ="Probability of Optimality (Preventing Infections)", x = "", fill = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p_res <- ggplot(prop_win_res, aes(y = Resistance, x = as.factor(Interventions))) + geom_bar(stat="identity") + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(y ="Probability of Optimality (Resistance Decrease)", x = "", fill = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggarrange(p_inf, p_res, ncol = 1, nrow = 2)

#Combination Plot 

p_comb <- ggplot(melt_combdata, aes(y = value, x = as.factor(Interventions), fill = variable)) + 
  geom_bar(stat="identity", position = position_dodge())  + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(y ="Probability of Optimality", x = "", fill = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

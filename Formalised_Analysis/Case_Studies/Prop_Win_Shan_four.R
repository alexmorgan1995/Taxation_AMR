library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Formalised_Analysis/Model_Output/Comparison_Sens")

# Import in Dataset -------------------------------------------------------

win_import <- readRDS("MDR_run_four.RDS")

for(i in seq_along(win_import)) {
  win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
}

# Altering Data Infections ------------------------------------------------

#Infections
#We want to minimise this - minimise the increase in total infections for each change in usage 
win_inf <- (win_import[,1:11])
#win_inf <- (win_import_pess[,1:10])
#win_inf <- (win_import_opt[,1:10])

win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("Infections" = colSums(win_inf_trans)/nrow(win_inf_trans),
                           "Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR1)","Single Tax (MR2)",
                                                         "Single Tax (LR)", 
                                                         "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                         "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                         "Diff Tax (5 Round)", "Diff Tax (6 Round)")))
prop_win_inf$Interventions <- factor(prop_win_inf$Interventions, levels = c(prop_win_inf$Interventions))

# Altering Data Resistance ------------------------------------------------

win_res <- (win_import[,12:22])

win_res_trans <- t(apply(win_res, 1, function(x) {
  val = max(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = colSums(win_res_trans)/nrow(win_res_trans),
                           "Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR1)","Single Tax (MR2)",
                                                         "Single Tax (LR)", 
                                                         "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                         "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                         "Diff Tax (5 Round)", "Diff Tax (6 Round)")))

prop_win_res$Interventions <- factor(prop_win_res$Interventions, levels = c(prop_win_res$Interventions))

# Altering Data Shannon's Index -------------------------------------------

win_shan <- round((win_import[,23:33]), 5)

win_shan[is.na(win_shan)] <- 0
win_shan <- win_shan[rowSums(win_shan[, -1]) > 0, ]

win_shan_trans <- t(apply(win_shan, 1, function(x) {
  val = max(x)
  x[is.na(x)] <- 0
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_shan <- data.frame("Shannon_Index" = colSums(win_shan_trans)/nrow(win_shan_trans),
                           "Interventions" = as.factor(c("Flat Tax", "Single Tax (HR)", "Single Tax (MR1)",
                                                         "Single Tax (MR2)","Single Tax (LR)", 
                                                         "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                                                         "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                                                         "Diff Tax (5 Round)", "Diff Tax (6 Round)")))

prop_win_shan$Interventions <- factor(prop_win_shan$Interventions, levels = c(prop_win_shan$Interventions))

# Combining All Together --------------------------------------------------

combdata <- prop_win_res; combdata$Infections <- prop_win_inf$Infections; combdata$Shannon <- prop_win_shan$Shannon_Index

melt_combdata <- melt(combdata, id.vars = "Interventions", measure.vars = c("Resistance", "Infections", "Shannon"))
melt_combdata$Interventions <- factor(melt_combdata$Interventions, levels = c(prop_win_res$Interventions))

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

p_shan <- ggplot(prop_win_shan, aes(y = Shannon_Index, x = as.factor(Interventions))) + geom_bar(stat="identity") + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(y ="Probability of Optimality (Highest SI)", x = "", fill = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggarrange(p_inf, p_res, p_shan, ncol = 1, nrow = 3)

#Combination Plot 

p_comb <- ggplot(melt_combdata, aes(y = value, x = as.factor(Interventions), fill = variable)) + 
  geom_bar(stat="identity", position = position_dodge())  + theme_bw() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  labs(y ="Probability of Intervention Winning", x = "", fill = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(values = c("red", "blue", "darkgreen"),labels = c("Resistance", "Infections", "Shannons Index"))

# Proportion of Wins HeatMap ----------------------------------------------

melt_combdata$value <- round(melt_combdata$value, digits = 3)

ggplot(melt_combdata, aes(Interventions, variable)) + theme_bw() +
  geom_tile(aes(fill = value)) + 
  facet_grid(variable ~ ., scales = "free_y") +
  geom_text(aes(label=value), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "Outcome Measure", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probabilty that \nIntervention Wins",
                               label.position = "bottom",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust=1))


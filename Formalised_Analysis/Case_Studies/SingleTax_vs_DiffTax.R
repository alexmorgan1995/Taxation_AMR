library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Formalised_Analysis/Model_Output/Comparison_Sens")

# Import in Dataset -------------------------------------------------------

win_import <- readRDS("MDR_run_v3_pess.RDS")

for(i in seq_along(win_import)) {
  win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
}

#Infections
win_inf <- (win_import[,1:10])
plot(density(win_inf[,1]))

#Testing for Normality
shapiro.test(win_inf[,1]) #Data is significantly different from a normal distribution - it is non normal
qqline(win_inf[,1], col = "steelblue", lwd = 2)

#Resistance
win_res <- (win_import[,11:20])
plot(density(win_res[,1]))

#Testing for Normality
shapiro.test(win_res[,1]) #Data is significantly different from a normal distribution - it is non-normal
qqline(win_res[,1], col = "steelblue", lwd = 2)

#Shannon's
win_shan <- round((win_import[,21:30]), 3)
win_shan[is.na(win_shan)] <- 0
win_shan <- win_shan[rowSums(win_shan[, -1]) > 0, ]

# Absolute Differences ----------------------------------------------------

#Infections

inc_inf <- win_inf[,c(2,5:10)]
inc_inf <- inc_inf[apply(inc_inf, 1, function(x) all(is.finite(x))), ]
keep <- Reduce(`&`, lapply(inc_inf, function(x) x >= quantile(x, .01) 
                           & x <= quantile(x, .99)))
inc_inf <- inc_inf[keep,]

m_inf <- melt(inc_inf, measure.vars = colnames(inc_inf))

#ggplot(m_inf, aes(x=variable, y=value, fill = variable)) + 
#  geom_violin() + theme_bw() + labs(y = "Relative Change from Single Tax (HR) Scenario", x = "") + 
#  scale_x_discrete(labels = c("Single Tax (HR)","Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) +
#  scale_fill_manual(values=c("darkgrey",mako(7)[-7])) + 
#  geom_boxplot(width=0.07, size = 0.4, color="white", outlier.shape = NA) + 
#  geom_signif(comparisons = list(c("single1_inf", "diff1_inf"),
#                                 c("single1_inf", "diff2_inf"),
#                                 c("single1_inf", "diff3_inf"),
#                                 c("single1_inf", "diff4_inf"),
#                                 c("single1_inf", "diff5_inf"),
#                                 c("single1_inf", "diff6_inf")),     
#              test = wilcox.test,
#              y_position = c(1.5, 1.6, 1.7, 1.8, 1.9, 2))

box_inf <- ggplot(m_inf, aes(x=variable, y=value, fill = variable)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Increase in Infections", x = "") + 
  scale_x_discrete(labels = c("Single Tax (HR)","Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) +
  scale_fill_manual(values=c("darkgrey", magma(8)[-c(1:2)])) + 
  geom_signif(comparisons = list(c("single1_inf", "diff1_inf"),
                                 c("single1_inf", "diff2_inf"),
                                 c("single1_inf", "diff3_inf"),
                                 c("single1_inf", "diff4_inf"),
                                 c("single1_inf", "diff5_inf"),
                                 c("single1_inf", "diff6_inf")),    
              map_signif_level=TRUE,    
              test = wilcox.test,
              
              tip_length = 0.02,
              vjust = 0.5,
              y_position = c(0.15, 0.175, 0.2, 0.225, 0.25, 0.275)) + coord_cartesian(ylim=c(-0.05, 0.5))

#Resistance

inc_res <- win_res[,c(2,5:10)]
inc_res <- inc_res[apply(inc_res, 1, function(x) all(is.finite(x))), ]
keep <- Reduce(`&`, lapply(inc_res, function(x) x >= quantile(x, .01) 
                           & x <= quantile(x, .99)))
inc_res <- inc_res[keep,]

m_res <- melt(inc_res, measure.vars = colnames(inc_res))

#ggplot(m_res, aes(x=variable, y=value, fill = variable)) + 
#  geom_violin() + theme_bw() + labs(y = "Relative Change in Intervention Efficacy (Resistance)", x = "") + 
#  scale_x_discrete(labels = c("Single Tax (HR)","Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) + 
#  scale_fill_manual(values=c("darkgrey",mako(7)[-7])) + 
#  geom_boxplot(width=0.07, size = 0.4, color="white", outlier.shape = NA) + 
#  geom_signif(comparisons = list(c("single1_inf", "diff1_inf"),
#                                 c("single1_inf", "diff2_inf"),
#                                 c("single1_inf", "diff3_inf"),
#                                 c("single1_inf", "diff4_inf"),
#                                 c("single1_inf", "diff5_inf"),
#                                 c("single1_inf", "diff6_inf")),   
#              map_signif_level=TRUE,
#              y_position = c(1.5, 1.6, 1.7, 1.8, 1.9, 2))

box_res <- ggplot(m_res, aes(x=variable, y=value, fill = variable)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Decrease in Resistance", x = "") + 
  scale_x_discrete(labels = c("Single Tax (HR)","Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) +
  scale_fill_manual(values=c("darkgrey",magma(8)[-c(1:2)])) + 
  geom_signif(comparisons = list(c("single1_res", "diff1_res"),
                                 c("single1_res", "diff2_res"),
                                 c("single1_res", "diff3_res"),
                                 c("single1_res", "diff4_res"),
                                 c("single1_res", "diff5_res"),
                                 c("single1_res", "diff6_res")),    
              map_signif_level=TRUE,    
              test = wilcox.test,
              tip_length = 0.02,
              vjust = 0.5,
              y_position = c(0.7, 0.775, 0.85, 0.925, 1, 1.75)) + coord_cartesian(ylim=c(-0.05, 1.5))

#Shannon's

diff_shan <- win_shan[,c(5:10)]
m_shan <- melt(diff_shan, measure.vars = colnames(diff_shan))

box_shan <-  ggplot(m_shan, aes(x=variable, y=value, fill = variable)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Shannon Index", x = "") + 
  scale_x_discrete(labels = c("Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) +
  scale_fill_manual(values=c(magma(8)[-c(1:2)])) 

#Kruskal Wallis

kruskal.test(value ~ variable, data = m_shan)

# Combine the Plots Together ----------------------------------------------

ggarrange(box_inf, box_res, labels= c("A", "B"), nrow = 2, ncol = 1)


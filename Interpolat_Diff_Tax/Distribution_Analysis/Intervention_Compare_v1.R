library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output")

# Import in Dataset -------------------------------------------------------

win_import <- readRDS("MDR_run_interpol.RDS")

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

#Average Antibiotics
win_avganti <- round((win_import[,31:40]), 3)


# Absolute Differences ----------------------------------------------------

#Infections

inc_inf <- win_inf[,c(1:10)]
inc_inf <- inc_inf[apply(inc_inf, 1, function(x) all(is.finite(x))), ]
keep <- Reduce(`&`, lapply(inc_inf, function(x) x >= quantile(x, .025) 
                           & x <= quantile(x, .975)))
inc_inf <- inc_inf[keep,]

m_inf <- melt(inc_inf, measure.vars = colnames(inc_inf))

ggplot(m_inf, aes(x=variable, y=value, fill = variable)) + 
  geom_violin() + theme_bw() + labs(y = "Relative Change from Single Tax (HR) Scenario", x = "") + 
  scale_x_discrete(labels = c("Flat Tax","Single Tax (HR)","Single Tax (MR)","Single Tax (LR)",
                              "Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", 
                              "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) +
  scale_fill_manual(values=c("darkgrey",mako(11)[-11])) + 
   stat_summary(fun.data=mean_sdl, size=0.5, 
                 geom="pointrange", color="red")

#Resistance

inc_res <- win_res[,c(1:10)]

inc_res <- inc_res[apply(inc_res, 1, function(x) all(is.finite(x))), ]
keep <- Reduce(`&`, lapply(inc_res, function(x) x >= quantile(x, .025) 
                           & x <= quantile(x, .975)))
inc_res <- inc_res[keep,]

m_res <- melt(inc_res, measure.vars = colnames(inc_res))

ggplot(m_res, aes(x=variable, y=value, fill = variable)) + 
  geom_violin() + theme_bw() + labs(y = "Relative Change from Single Tax (HR) Scenario", x = "") + 
  scale_x_discrete(labels = c("Flat Tax","Single Tax (HR)","Single Tax (MR)","Single Tax (LR)",
                              "Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", 
                              "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) +
  scale_fill_manual(values=c("darkgrey",mako(11)[-11])) + 
  stat_summary(fun.data=mean_sdl, size=0.5, 
               geom="pointrange", color="red")

#Average Antibiotic Available 

inc_avganti <- win_avganti[,c(1:10)]
m_avganti <- melt(inc_avganti, measure.vars = colnames(inc_avganti))

ggplot(m_avganti, aes(x=variable, y=value, fill = variable)) + 
  geom_violin() + theme_bw() + labs(y = "Relative Change from Single Tax (HR) Scenario", x = "") + 
  scale_x_discrete(labels = c("Flat Tax","Single Tax (HR)","Single Tax (MR)","Single Tax (LR)",
                              "Diff Tax (1 Rd)", "Diff Tax (2 Rds)", "Diff Tax (3 Rds)", 
                              "Diff Tax (4 Rds)", "Diff Tax (5 Rds)", "Diff Tax (6 Rds)")) +
  scale_fill_manual(values=c("darkgrey",mako(11)[-11])) + 
  stat_summary(fun.data=mean_sdl, size=0.5, 
               geom="pointrange", color="red")




#Shannon's

diff_shan <- win_shan[,c(1,5:10)]
m_shan <- melt(diff_shan, measure.vars = colnames(diff_shan))

box_shan <- ggplot(m_shan, aes(x=variable, y=value, fill = variable)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Shannon Index", x = "") 

kruskal.test(value ~ variable, data = m_shan)

# Combine the Plots Together ----------------------------------------------

ggarrange(box_inf, box_res, 
          box_avganti, labels= c("A", "B", "C"), nrow = 3, ncol = 1)




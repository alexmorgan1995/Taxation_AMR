library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/New")

# Import in Dataset -------------------------------------------------------

win_import <- readRDS("MDR_run_interpol_v1.RDS")

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

box_inf <- ggplot(m_inf, aes(x=variable, y=value, fill = variable)) + coord_cartesian(ylim=c(-0.125, 0.175)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Increase in Infections", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=11), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_discrete(labels= c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                             "Single Tax (LR)", 
                             "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                             "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                             "Diff Tax (5 Round)", "Diff Tax (6 Round)"))

pairwise.wilcox.test(m_inf$value, m_inf$variable,
                     p.adjust.method = "BH")

#Resistance

inc_res <- win_res[,c(1:10)]
m_res <- melt(inc_res, measure.vars = colnames(inc_res))

box_res <- ggplot(m_res, aes(x=variable, y=value, fill = variable)) + coord_cartesian(ylim=c(-0.5, 1.25)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Decrease in Resistance", x = "")  + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=11), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_discrete(labels= c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                             "Single Tax (LR)", 
                             "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                             "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                             "Diff Tax (5 Round)", "Diff Tax (6 Round)"))

kruskal.test(value ~ variable, data = m_inf)

#Average Antibiotic Available 

inc_avganti <- win_avganti[,c(1:10)]
m_avganti <- melt(inc_avganti, measure.vars = colnames(inc_avganti))

box_avganti <-  ggplot(m_avganti, aes(x=variable, y=value, fill = variable)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Average Antibiotic Available", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=11), axis.title.x= element_text(size=11), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_discrete(labels= c("Flat Tax", "Single Tax (HR)", "Single Tax (MR)",
                             "Single Tax (LR)", 
                             "Diff Tax (1 Round)", "Diff Tax (2 Round)",
                             "Diff Tax (3 Round)", "Diff Tax (4 Round)", 
                             "Diff Tax (5 Round)", "Diff Tax (6 Round)"))

kruskal.test(value ~ variable, data = m_avganti)



#Shannon's

diff_shan <- win_shan[,c(1,5:10)]
m_shan <- melt(diff_shan, measure.vars = colnames(diff_shan))

box_shan <- ggplot(m_shan, aes(x=variable, y=value, fill = variable)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Shannon Index", x = "") 

kruskal.test(value ~ variable, data = m_shan)

# Combine the Plots Together ----------------------------------------------

ggarrange(box_res, box_inf, 
          box_avganti, labels= c("A", "B", "C"), nrow = 3, ncol = 1)




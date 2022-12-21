library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
rm(list=ls())

# Import in Dataset -------------------------------------------------------

tax_data <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/taxlist.RDS")

# Resistance --------------------------------------------------------------

tax_dist <- tax_data
m_tax_dist <- melt(tax_dist, measure.vars = colnames(tax_dist))

test_stat <- pairwise.wilcox.test(m_tax_dist$value, m_tax_dist$variable,
                     p.adjust.method = "bonferroni")

#Box Plot
box_tax <- ggplot(m_tax_dist, aes(x=variable, y=(value/1000000000)/200, fill = variable)) + coord_cartesian(ylim=c(0, 3.5)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Average Yearly Revenue ($ 10 Billion)", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), 
        axis.text=element_text(size=11), 
        axis.title =element_text(size=12),  plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)")) 

# Output -----------------------------------------------------------------

ggsave(box_tax, filename = "tax_distribution.png", dpi = 300, width = 9, height = 4, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

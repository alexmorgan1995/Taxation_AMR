library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
rm(list=ls())

# Import in Dataset -------------------------------------------------------

tax_data <- readRDS("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Money_Analysis/taxlist.RDS")

# Resistance --------------------------------------------------------------

tax_dist <- tax_data
m_tax_dist <- melt(tax_dist, measure.vars = colnames(tax_dist))
m_tax_dist$country <- rep(c("USA", "China"), each = 1000, times = 10)
  
test_stat <- pairwise.wilcox.test(m_tax_dist$value, m_tax_dist$variable,
                     p.adjust.method = "bonferroni")

# Faceted Box Plot --------------------------------------------------------

#Need this function to be able to limit the axis when faceting (due to the outliers)
calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

#Box Plot
box_tax <- ggplot(m_tax_dist, aes(x=variable, y=(value/1000000000)/200, fill=variable)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + theme_bw() +
  facet_wrap(~country, scales="free", ncol = 1 , nrow = 2)  + labs(y = "Average Yearly Revenue ($ 10 Billion)", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), 
        axis.text=element_text(size=11), 
        axis.title =element_text(size=12),  plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels= c("FT", "ST (HR)", 
                             "ST (MR)", "ST (LR)", 
                             "DT (1Rd)",
                             "DT (2Rd)", 
                             "DT (3Rd)", 
                             "DT (4Rd)", 
                             "DT (5Rd)",
                             "DT (6Rd)")) + guides(fill = "none") 

# Average Box Plot --------------------------------------------------------

#Box Plot
country_box_tax <- ggplot(m_tax_dist, aes(x=country, y=(value/1000000000)/200, fill = country)) + coord_cartesian(ylim=c(0, 3)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Average Yearly Revenue ($ 10 Billion)", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), 
        axis.text=element_text(size=11), 
        axis.title =element_text(size=12),  plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) + 
  scale_y_continuous(expand = c(0,0)) 

# Output -----------------------------------------------------------------

tax_plot <- ggarrange(box_tax, country_box_tax, nrow = 1, ncol = 2, widths = c(1, 0.4), labels = c("A", "B"), 
          font.label = list(size = 18, color = "black", face = "bold"))

ggplot_build(country_box_tax)$data[[1]]

ggsave(tax_plot, filename = "tax_distribution.png", dpi = 300, width = 12, height = 6, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

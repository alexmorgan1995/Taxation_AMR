library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Extra_Anal")

# Import in the RDS -------------------------------------------------------

#So the way in which this is built, is that we have a 1000 sets of uncertainty analysis
#With 100 taxation rates
#In each taxation rate there is an uncertainty analysis 

threshold_analysis <- readRDS("uncert_thresh.RDS")
thresh_dataframe_tax <- data.frame(matrix(ncol = 100, nrow = 1000))
thresh_dataframe_bans <- data.frame(matrix(ncol = 100, nrow = 1000))

for(i in 1:100) {
  test <- threshold_analysis[[i]]
  thresh_dataframe_tax[,i] <- unlist(sapply(test,function(x) x[3]))
  thresh_dataframe_bans[,i] <- unlist(sapply(test,function(x) x[4]))
}

colnames(thresh_dataframe_tax) <- seq(1,100)
colnames(thresh_dataframe_bans) <- seq(1,100)

thresh_dataframe_tax[thresh_dataframe_tax == -1000] = NA
thresh_dataframe_bans[thresh_dataframe_bans == -1000] = NA

thresh_dataframe_tax <- thresh_dataframe_tax[,c(1, seq(5,100, by = 5))]
thresh_dataframe_bans <- thresh_dataframe_bans[,c(1, seq(5,100, by = 5))]

#Taxation
melt_tax <- melt(thresh_dataframe_tax)

taxation <- ggplot(melt_tax, aes(x = variable, y = value, fill  = value))+
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "grey") + coord_cartesian(ylim=c(-1, 0.3)) + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title=element_text(size=11), axis.title.x=element_text(size=11), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(size=11),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) + 
  scale_x_discrete(labels= c(1,seq(5,100, by = 5)), name = "Taxation (%)") + 
  scale_y_continuous(name = "Change in Resistance under Curtailment") + 
  geom_hline(yintercept = median(thresh_dataframe_tax[,1], na.rm = T))

median(thresh_dataframe_tax[,13], na.rm = T)

ggplot_build(taxation)$data

#Bans
melt_bans <- melt(thresh_dataframe_bans)

bans <- ggplot(melt_bans, aes(x = variable, y = value, fill  = value))+
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "grey") + coord_cartesian(ylim=c(-1.1, 0.5)) + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title=element_text(size=11), axis.title.x=element_text(size=11), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(size=11), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) + 
  scale_x_discrete(labels= c(1,seq(5,100, by = 5)), name = "Taxation (%)") + 
  scale_y_continuous(name = "Change in Resistance under Curtailment") + 
  geom_hline(yintercept = median(thresh_dataframe_bans[,1], na.rm = T))

ggplot_build(bans)$data

ggsave(taxation, filename = "thresh_tax.png", dpi = 300, width = 11, height = 4, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

ggsave(bans, filename = "thresh_bans.png", dpi = 300, width = 11, height = 4, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")



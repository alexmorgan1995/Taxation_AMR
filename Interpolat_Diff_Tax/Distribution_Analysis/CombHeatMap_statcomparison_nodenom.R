library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
rm(list=ls())
setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/Bans/")

# Import in Dataset -------------------------------------------------------

win_import_change <- readRDS("MDR_run_Base_tax_nodenom.RDS"); win_import <- win_import_change

win_import[win_import == -1000] <- NA

for(i in seq_along(win_import)) {
  win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))
}

#Infections
win_inf <- (win_import[,1:13])
#plot(density(win_inf[,1]))

#Testing for Normality
shapiro.test(win_inf[,1]) #Data is significantly different from a normal distribution - it is non normal
#qqline(win_inf[,1], col = "steelblue", lwd = 2)

#Resistance
win_res <- (win_import[,14:26])

#Testing for Normality
shapiro.test(win_res[,1]) #Data is significantly different from a normal distribution - it is non-normal
#qqline(win_res[,1], col = "steelblue", lwd = 2)

#Shannon's
win_shan <- round((win_import[,27:39]), 3)
win_shan[is.na(win_shan)] <- 0
win_shan <- win_shan[rowSums(win_shan[, -1]) > 0, ]

#Average Antibiotics
win_avganti <- round((win_import[,40:52]), 3)

# Resistance --------------------------------------------------------------

inc_res <- win_res
m_res <- melt(inc_res, measure.vars = colnames(inc_res))
m_res$factors <- c(rep("Taxation", 10000), rep("Ban", 3000))
m_res$factors <- factor(m_res$factors, levels = c("Taxation", "Ban"))

test_stat_res <- pairwise.wilcox.test(m_res$value, m_res$variable,
                     p.adjust.method = "bonferroni")

#Intervention Failure
prop_vec <- data.frame("intervention" = unique(m_res$variable),
                       "prop_inc" = NA,
                       "95_quant" = NA)

m_win_import_change <- melt(win_import_change, measure.vars = colnames(win_import_change))

for(i in 1:13) {
  prop_1000 <- win_import_change[,c(14:26)[i]]
  prop_vec[i,2] <- length(prop_1000[prop_1000 == -1000])/sum(!is.na(prop_1000))
}

p <- ggplot(m_res, aes(x=factor(variable), y=value)) + geom_boxplot() 
prop_vec[,3] <- ggplot_build(p)$data[[1]]$ymax

#Probability of Winning
win_res_trans <- t(apply(inc_res, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = (colSums(win_res_trans, na.rm = T)/nrow(win_res_trans))*100,
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (MR)", "Ban (LR)")))

#Tranformed Win
prop_win_res$Resistance <- prop_win_res$Resistance*(1-prop_vec$prop_inc)
prop_win_res$Resistance <- round(prop_win_res$Resistance/sum(prop_win_res$Resistance)*100,1)

prop_win_res$Interventions <- factor(prop_win_res$Interventions, levels = c(prop_win_res$Interventions))
prop_win_res$factors <- c(rep("Taxation", 10) , rep("Ban", 3))
prop_win_res$factors <- factor(prop_win_res$factors, levels = c("Taxation", "Ban"))
prop_win_res$Color <- "black" 
prop_win_res[prop_win_res$Resistance == max(prop_win_res$Resistance),4] <- "white"

#Win Heat Map
win_res_p <- ggplot(prop_win_res, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Resistance)) + 
  facet_grid(. ~ factors, scales = "free") +
  geom_text(aes(label=Resistance, color = Interventions)) + 
  scale_colour_manual(values=prop_win_res$Color) +
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins (%)",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1), color = "none") + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.75), "cm"))

#Box Plot
box_res <- ggplot(m_res, aes(x=variable, y=value/(10300-3000), fill = variable, alpha = variable)) + coord_cartesian(ylim=c(-0.35, 0.1)) +
  facet_grid(. ~ factors, scales = "free", space = "free") +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "red") + theme_bw() + labs(y = "Change in Resistance under Curtailment", x = "") + 
  scale_alpha_manual(values=  prop_vec$prop_inc) +
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title=element_text(size=11), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)",
                             "Ban (HR)", "Ban (MR)", "Ban (LR)")) 

ggplot_build(box_res)$data

comb_res <- ggarrange(box_res, win_res_p, ncol =1, nrow= 2, heights = c(1, 0.5), align = "v")

# Absolute Differences ----------------------------------------------------

#Infections Distribution

inc_inf <- win_inf
#inc_inf <- inc_inf[apply(inc_inf, 1, function(x) all(is.finite(x))), ]
#keep <- Reduce(`&`, lapply(inc_inf, function(x) x >= quantile(x, .025) 
#                           & x <= quantile(x, .975)))
#inc_inf <- inc_inf[keep,]

m_inf <- melt(inc_inf, measure.vars = colnames(inc_inf))
m_inf$factors <- c(rep("Taxation", 10000), rep("Ban", 3000))
m_inf$factors <- factor(m_inf$factors, levels = c("Taxation", "Ban"))

test_stat_inf <- pairwise.wilcox.test(m_inf$value, m_inf$variable,
                                  p.adjust.method = "bonferroni")

#Infections Wins
inc_inf[is.na(inc_inf)] <- 0
win_inf_trans <- t(apply(inc_inf, 1, function(x) {
  val = min(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("Infections" = (round(colSums(win_inf_trans)/nrow(win_inf_trans),3))*100,
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)",
                                                         "Ban (HR)", "Ban (MR)", "Ban (LR)")))

#Tranformed Win
prop_win_inf$Infections <- prop_win_inf$Infections*(1-prop_vec$prop_inc)
prop_win_inf$Infections <- round(prop_win_inf$Infections/sum(prop_win_inf$Infections)*100,1)

prop_win_inf$Interventions <- factor(prop_win_inf$Interventions, levels = c(prop_win_inf$Interventions))
prop_win_inf$factors <- c(rep("Taxation", 10), rep("Ban", 3))
prop_win_inf$factors <- factor(prop_win_inf$factors, levels = c("Taxation", "Ban"))
prop_win_inf$Color <- "black" 
prop_win_inf[prop_win_inf$Infections == max(prop_win_inf$Infections),4] <- "white"

#Win Heat Map
win_inf_p <- ggplot(prop_win_inf, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Infections)) + 
  facet_grid(. ~ factors, scales = "free", space = "free") +
  geom_text(aes(label=Infections, color = Interventions)) + 
  scale_colour_manual(values=prop_win_inf$Color) +
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins (%)",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1), color = "none") + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.75), "cm"))

#Box Plot
box_inf <- ggplot(m_inf, aes(x=variable, y=value/(10300-3000), fill = variable, alpha = variable)) + coord_cartesian(ylim=c(-0.05, 0.075)) + 
  facet_grid(. ~ factors, scales = "free", space = "free") +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "red") + theme_bw() + labs(y = "Change in Infections under Curtailment", x = "") + 
  scale_alpha_manual(values=  prop_vec$prop_inc) +
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title=element_text(size=11), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)",
                             "Ban (HR)", "Ban (MR)", "Ban (LR)"))

ggplot_build(box_inf)$data

comb_inf <- ggarrange(box_inf, win_inf_p, ncol =1, nrow= 2, heights = c(1, 0.5), align = "v")

# Average Antibiontics Available ------------------------------------------

inc_avganti <- win_avganti
m_avganti <- melt(inc_avganti, measure.vars = colnames(inc_avganti))
m_avganti$factors <- c(rep("Taxation", 10000), rep("Ban", 3000))
m_avganti$factors <- factor(m_avganti$factors, levels = c("Taxation", "Ban"))


test_stat_avganti <- pairwise.wilcox.test(m_avganti$value, m_avganti$variable,
                                  p.adjust.method = "bonferroni")

inc_avganti[is.na(inc_avganti)] <- 0
inc_avganti <- inc_avganti[rowSums(inc_avganti[, -1]) > 0, ]

win_avganti_trans <- t(apply(inc_avganti, 1, function(x) {
  val = max(x)
  x[is.na(x)] <- 0
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_avganti <- data.frame("Average_Anti" = (round(colSums(win_avganti_trans)/nrow(win_avganti_trans),3))*100,
                               "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                               "ST (LR)", 
                                                               "DT (1Rd)", "DT (2Rd)",
                                                               "DT (3Rd)", "DT (4Rd)", 
                                                               "DT (5Rd)", "DT (6Rd)",
                                                               "Ban (HR)", "Ban (MR)", "Ban (LR)")))

#Tranformed Win
prop_win_avganti$Average_Anti <- prop_win_avganti$Average_Anti*(1-prop_vec$prop_inc)
prop_win_avganti$Average_Anti <- round(prop_win_avganti$Average_Anti/sum(prop_win_avganti$Average_Anti)*100,1)

prop_win_avganti$Interventions <- factor(prop_win_avganti$Interventions, levels = c(prop_win_avganti$Interventions))
prop_win_avganti$factors <- c(rep("Taxation", 10), rep("Ban", 3))
prop_win_avganti$factors <- factor(prop_win_avganti$factors, levels = c("Taxation", "Ban"))
prop_win_avganti$Color <- "black" 
prop_win_avganti[prop_win_avganti$Average_Anti == max(prop_win_avganti$Average_Anti),4] <- "white"

#Win Heat Map

win_avganti_p <- ggplot(prop_win_avganti, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Average_Anti)) + 
  facet_grid(. ~ factors, scales = "free", space = "free") +
  geom_text(aes(label=Average_Anti, color = Interventions)) + 
  scale_colour_manual(values=prop_win_avganti$Color) +
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins (%)",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1), color = "none") + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.75), "cm"))

#Box Plot

box_avganti <- ggplot(m_avganti, aes(x=variable, y=value, fill = variable, alpha = variable)) + coord_cartesian(ylim=c(-0.0001, 3)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "red") +
  facet_grid(. ~ factors, scales = "free_x", space = "free")  + theme_bw() + labs(y = "Number of Available Antibiotics", x = "")  + 
  scale_alpha_manual(values =  prop_vec$prop_inc) +
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title=element_text(size=11), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,2), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 11, colour = "black", face="bold")) 

comb_avg_anti <- ggarrange(box_avganti, win_avganti_p, ncol = 1, nrow= 2, heights = c(1, 0.5), align = "v")

# Combine the Plots Together ----------------------------------------------

test <- ggarrange(comb_res, comb_inf, 
          comb_avg_anti, labels= c("A", "B", "C"), font.label=list(color="black",size=20) ,nrow = 3, ncol = 1, align="v",
          heights = c(0.1, 0.1, 0.1), common.legend = T)

ggsave(test, filename = "run_nodenom.png", dpi = 300, width = 11, height = 14, units = "in",
       path = "/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Figures/")

box_avganti <- box_avganti + theme(axis.text.x = element_text(size=12)) 

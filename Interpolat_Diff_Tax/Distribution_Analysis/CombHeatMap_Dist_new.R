library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")
rm(list = ls())
setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/New_v2")

# Import in Dataset -------------------------------------------------------

win_import_change <- readRDS("MDR_run_interpol_new.RDS"); win_import <- win_import_change
#win_import_change <- readRDS("MDR_run_interpol_biasPED.RDS"); win_import <- win_import_change
#win_import_change <- readRDS("MDR_run_interpol_realPED.RDS"); win_import <- win_import_change
#win_import_change <- readRDS("MDR_run_interpol_25.RDS"); win_import <- win_import_change
#win_import_change <- readRDS("MDR_run_interpol_75.RDS"); win_import <- win_import_change

win_import[win_import == -1000] <- NA

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

#Testing for Normality
shapiro.test(win_res[,1]) #Data is significantly different from a normal distribution - it is non-normal
qqline(win_res[,1], col = "steelblue", lwd = 2)

#Shannon's
win_shan <- round((win_import[,21:30]), 3)
win_shan[is.na(win_shan)] <- 0
win_shan <- win_shan[rowSums(win_shan[, -1]) > 0, ]

#Average Antibiotics
win_avganti <- round((win_import[,31:40]), 3)


# Resistance --------------------------------------------------------------

win_res <- (win_import[,11:20])
inc_res <- win_res[,c(1:10)]

m_res <- melt(inc_res, measure.vars = colnames(inc_res))

prop_vec <- data.frame("intervention" = unique(m_res$variable),
                       "prop_inc" = NA,
                       "95_quant" = NA )

m_win_import_change <- melt(win_import_change, measure.vars = colnames(win_import_change))

for(i in 1:10) {
  prop_1000 <- win_import_change[,i+10]
  prop_vec[i,2] <- length(prop_1000[prop_1000 == -1000])/sum(!is.na(prop_1000))
}

p <- ggplot(m_res, aes(x=factor(variable), y=value)) + geom_boxplot() 
prop_vec[,3] <- ggplot_build(p)$data[[1]]$ymax

#win_res[is.na(win_res)] <- 0
win_res_trans <- t(apply(win_res, 1, function(x) {
  val = max(x, na.rm = T)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = colSums(win_res_trans, na.rm = T)/nrow(win_res_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)")))

prop_win_res$Interventions <- factor(prop_win_res$Interventions, levels = c(prop_win_res$Interventions))
prop_win_res$Color <- "black" 
prop_win_res[prop_win_res$Resistance == max(prop_win_res$Resistance),3] <- "white"

#Win Heat Map

win_res_p <- ggplot(prop_win_res, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Resistance)) + 
  geom_text(aes(label=Resistance, color = Interventions)) + 
  scale_colour_manual(values=prop_win_res$Color) +
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
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
box_res <- ggplot(m_res, aes(x=variable, y=value, fill = variable, alpha = variable)) + coord_cartesian(ylim=c(-1, 6)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "red") + theme_bw() + labs(y = "Decrease in Resistance (%)", x = "") + 
  scale_alpha_manual(values=  prop_vec$prop_inc) + guides(fill="none") +
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)")) +
  geom_signif(comparisons = list(c("singleMR_res", "flat_res"),
                                 c("singleMR_res", "singleHR_res"),
                                 c("singleMR_res", "singleLR_res"),
                                 c("singleMR_res", "diff1_res"),
                                 c("singleMR_res", "diff2_res"),
                                 c("singleMR_res", "diff3_res"),
                                 c("singleMR_res", "diff4_res"),
                                 c("singleMR_res", "diff5_res"),
                                 c("singleMR_res", "diff6_res")), 
              map_signif_level=TRUE, show.legend = F,
              tip_length = 0,
              y_position = c(-32.5, -32, -31.5, -31, -30.5, -30, -29.5, -29, -28.5))

comb_res <- ggarrange(box_res, win_res_p, ncol =1, nrow= 2, heights = c(1, 0.6), align = "v")

# Absolute Differences ----------------------------------------------------

#Infections Distribution

inc_inf <- win_inf[,c(1:10)]

m_inf <- melt(inc_inf, measure.vars = colnames(inc_inf))

#Infections Wins

win_inf[is.na(win_inf)] <- 0
win_inf_trans <- t(apply(win_inf, 1, function(x) {
  val = min(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_inf <- data.frame("Infections" = colSums(win_inf_trans)/nrow(win_inf_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)")))

prop_win_inf$Interventions <- factor(prop_win_inf$Interventions, levels = c(prop_win_inf$Interventions))
prop_win_inf$Color <- "black" 
prop_win_inf[prop_win_inf$Infections == max(prop_win_inf$Infections),3] <- "white"

#Win Heat Map
win_inf_p <- ggplot(prop_win_inf, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Infections)) + 
  geom_text(aes(label=Infections, color = Interventions)) + 
  scale_colour_manual(values=prop_win_inf$Color) +
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
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
box_inf <- ggplot(m_inf, aes(x=variable, y=value, fill = variable, alpha = variable)) +  coord_cartesian(ylim=c(-0.125, 4)) + 
  scale_alpha_manual(values=  prop_vec$prop_inc) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "red") + theme_bw() + labs(y = "Increase in Infections (%)", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)")) +
  geom_signif(comparisons = list(c("singleLR_inf", "flat_inf"),
                                 c("singleLR_inf", "singleHR_inf"),
                                 c("singleLR_inf", "singleMR_inf"),
                                 c("singleLR_inf", "diff1_inf"),
                                 c("singleLR_inf", "diff2_inf"),
                                 c("singleLR_inf", "diff3_inf"),
                                 c("singleLR_inf", "diff4_inf"),
                                 c("singleLR_inf", "diff5_inf"),
                                 c("singleLR_inf", "diff6_inf")), 
              map_signif_level=TRUE, show.legend = F,
              tip_length = 0,
              y_position = c(-33.2, -32.8, -32.4, -32, -31.6, -31.2, -30.8, -30.4, -30))

comb_inf <- ggarrange(box_inf, win_inf_p, ncol =1, nrow= 2, heights = c(1, 0.6), align = "v")

# Average Antibiontics Available ------------------------------------------

inc_avganti <- win_avganti[,c(1:10)]
m_avganti <- melt(inc_avganti, measure.vars = colnames(inc_avganti))

win_avganti[is.na(win_avganti)] <- 0
win_avganti <- win_avganti[rowSums(win_avganti[, -1]) > 0, ]

win_avganti_trans <- t(apply(win_avganti, 1, function(x) {
  val = max(x)
  x[is.na(x)] <- 0
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_avganti <- data.frame("Average_Anti" = round(colSums(win_avganti_trans)/nrow(win_avganti_trans),3),
                               "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                             "ST (LR)", 
                                                             "DT (1Rd)", "DT (2Rd)",
                                                             "DT (3Rd)", "DT (4Rd)", 
                                                             "DT (5Rd)", "DT (6Rd)")))

prop_win_avganti$Interventions <- factor(prop_win_avganti$Interventions, levels = c(prop_win_avganti$Interventions))
prop_win_avganti$Color <- "black" 
prop_win_avganti[prop_win_avganti$Average_Anti == max(prop_win_avganti$Average_Anti),3] <- "white"

#Win Heat Map
win_avganti_p <- ggplot(prop_win_avganti, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Average_Anti)) + 
  geom_text(aes(label=Average_Anti, color = Interventions)) + 
  scale_colour_manual(values=prop_win_avganti$Color) + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1), color = "none") + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.35), "cm"))

#Box Plot
box_avganti <- ggplot(m_avganti, aes(x=variable, y=value, fill = variable, alpha = variable)) + coord_cartesian(ylim=c(-0.001, 3)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "red") + theme_bw() + labs(y = "Average Antibiotics Available", x = "") + 
  scale_alpha_manual(values=  prop_vec$prop_inc) +
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,2), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)"))

comb_avg_anti <- ggarrange(box_avganti, win_avganti_p, ncol =1, nrow= 2, heights = c(1, 0.6), align = "v")

# Shannon's ---------------------------------------------------------------

diff_shan <- win_shan[,c(1:10)]
m_shan <- melt(diff_shan, measure.vars = colnames(diff_shan))

win_shan[is.na(win_shan)] <- 0
win_shan <- win_shan[rowSums(win_shan[, -1]) > 0, ]

win_shan_trans <- t(apply(win_shan, 1, function(x) {
  val = max(x)
  x[is.na(x)] <- 0
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_shan <- data.frame("Shannon_Index" = round(colSums(win_shan_trans)/nrow(win_shan_trans),3),
                            "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                          "ST (LR)", 
                                                          "DT (1Rd)", "DT (2Rd)",
                                                          "DT (3Rd)", "DT (4Rd)", 
                                                          "DT (5Rd)", "DT (6Rd)")))

prop_win_shan$Interventions <- factor(prop_win_shan$Interventions, levels = c(prop_win_shan$Interventions))
prop_win_shan$Color <- "black" 
prop_win_shan[prop_win_shan$Shannon_Index == max(prop_win_shan$Shannon_Index),3] <- "white"

#Win Heat Map
win_shan_p <- ggplot(prop_win_shan, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Shannon_Index)) + 
  geom_text(aes(label=Shannon_Index, color = Interventions)) + 
  scale_colour_manual(values=prop_win_shan$Color) + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probability that Intervention Wins", 
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1), color = "none") + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.8), "cm"))

#Box Plot
box_shan <- ggplot(m_shan, aes(x=variable, y=value, fill = variable, alpha = variable)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, fill = "red") + theme_bw() + labs(y = "Shannon's Index", x = "") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_alpha_manual(values=  prop_vec$prop_inc) +
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)"))


comb_shan <- ggarrange(box_shan, win_shan_p, ncol =1, nrow= 2, heights = c(1, 0.6), align = "v")

# Combine the Plots Together ----------------------------------------------

test <- ggarrange(comb_res, comb_inf,
                  comb_avg_anti, labels= c("A", "B", "C"), font.label=list(color="black",size=20) ,nrow = 3, ncol = 1, align="hv",
                  heights = c(0.1, 0.1, 0.1), common.legend = T)


ggsave(test, filename = "test_full.png", dpi = 300, width = 9, height = 12, units = "in",
       path = "/Users/amorgan/Desktop")

ggsave(comb_shan, filename = "shan_full.png", dpi = 300, width = 9, height = 3.5, units = "in",
       path = "/Users/amorgan/Desktop")


setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Interpolat_Diff_Tax/Euler_Run/Model_Output/New_V1")

# Import in Dataset -------------------------------------------------------

win_import <- readRDS("MDR_run_interpol_avgresinf.RDS")

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

#Infections Distribution

inc_inf <- win_inf[,c(1:10)]
inc_inf <- inc_inf[apply(inc_inf, 1, function(x) all(is.finite(x))), ]
keep <- Reduce(`&`, lapply(inc_inf, function(x) x >= quantile(x, .025) 
                           & x <= quantile(x, .975)))
inc_inf <- inc_inf[keep,]
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

#Win Heat Map
win_inf_p <- ggplot(prop_win_inf, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Infections)) + 
  geom_text(aes(label=Infections), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probabilty that Intervention Wins",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      plot.margin = unit(c(0,1,0,1.75), "cm"))

#Box Plot
box_inf <- ggplot(m_inf, aes(x=variable, y=value, fill = variable)) + coord_cartesian(ylim=c(-0.125, 0.175)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Increase in Infections (%)", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)"))

comb_inf <- ggarrange(box_inf, win_inf_p, ncol =1, nrow= 2, heights = c(1, 0.6), align = "v")


# Resistance --------------------------------------------------------------

inc_res <- win_res[,c(1:10)]
m_res <- melt(inc_res, measure.vars = colnames(inc_res))

win_res <- (win_import[,11:20])
win_res[is.na(win_res)] <- 0
win_res_trans <- t(apply(win_res, 1, function(x) {
  val = max(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

prop_win_res <- data.frame("Resistance" = colSums(win_res_trans)/nrow(win_res_trans),
                           "Interventions" = as.factor(c("FT", "ST (HR)", "ST (MR)",
                                                         "ST (LR)", 
                                                         "DT (1Rd)", "DT (2Rd)",
                                                         "DT (3Rd)", "DT (4Rd)", 
                                                         "DT (5Rd)", "DT (6Rd)")))

prop_win_res$Interventions <- factor(prop_win_res$Interventions, levels = c(prop_win_res$Interventions))

#Win Heat Map
win_res_p <- ggplot(prop_win_res, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Resistance)) + 
  geom_text(aes(label=Resistance), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probabilty that Intervention Wins",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.75), "cm"))

#Box Plot
box_res <- ggplot(m_res, aes(x=variable, y=value, fill = variable)) + coord_cartesian(ylim=c(-0.5, 1.25)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Decrease in Resistance (%)", x = "") + 
  theme(legend.position= "bottom", legend.text=element_text(size=11), legend.title =element_text(size=12), axis.text=element_text(size=11), 
        axis.title.y=element_text(size=12), axis.title.x= element_blank(), plot.margin = unit(c(0.3,1,0,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  scale_x_discrete(labels= c("FT", "ST (HR)", "ST (MR)",
                             "ST (LR)", 
                             "DT (1Rd)", "DT (2Rd)",
                             "DT (3Rd)", "DT (4Rd)", 
                             "DT (5Rd)", "DT (6Rd)"))

comb_res <- ggarrange(box_res, win_res_p, ncol =1, nrow= 2, heights = c(1, 0.6), align = "v")

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

#Win Heat Map
win_avganti_p <- ggplot(prop_win_avganti, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Average_Anti)) + 
  geom_text(aes(label = Average_Anti), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probabilty that Intervention Wins",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.35), "cm"))

#Box Plot
box_avganti <- ggplot(m_avganti, aes(x=variable, y=value, fill = variable)) + coord_cartesian(ylim=c(-0.001, 3)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Average Antibiotics Available", x = "") + 
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

#Win Heat Map
win_shan_p <- ggplot(prop_win_shan, aes(Interventions, "")) + theme_bw() +
  geom_tile(aes(fill = Shannon_Index)) + 
  geom_text(aes(label = Shannon_Index), color = "black") + 
  scale_fill_distiller(palette ="Blues", direction = 1) +
  scale_x_discrete(name = "", expand = c(0, 0))  +   
  scale_y_discrete(name = "", expand = c(0, 0)) + 
  guides(fill = guide_colorbar(title = "Probabilty that Intervention Wins",
                               title.position = "left", title.vjust = 1,
                               # draw border around the legend
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1)) + 
  theme(strip.background = element_blank(), axis.text=element_text(size=11),
        strip.text = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 0, hjust=0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,1,0,1.8), "cm"))

#Box Plot
box_shan <- ggplot(m_shan, aes(x=variable, y=value, fill = variable)) + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() + labs(y = "Shannon's Index", x = "") + 
  scale_y_continuous(expand = c(0,0)) +
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



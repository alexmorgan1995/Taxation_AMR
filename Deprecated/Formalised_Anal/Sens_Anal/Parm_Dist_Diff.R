library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/Theoretical_Analysis/Formalised_Anal/Sens_Anal/PED_Scen/")


# Parameter Import --------------------------------------------------------

win_import <- readRDS("comb_data.RDS")
for(i in seq_along(win_import)) win_import[[i]] <- as(win_import[[i]], class(win_import[[i]][[1]]))


# Parm Set ----------------------------------------------------------------

low_parm <- c(0, #beta
              1/50, #mu_w
              1/50, #mu_r
              1/50, #mu_t
              0, #eta_wr
              0, #eta_wr
              0.5, #c1
              0.5, #c2
              0.5, #c3
              0, #rho
              0) #baseline tax

high_parm <- c(10, #beta
               1/5, #mu_w
               1/5, #mu_r
               1/5, #mu_t
               1, #eta_wr
               1, #eta_wr
               1, #c1
               1, #c2
               1, #c3
               1, #rho
               1) #baseline tax

# Alter Dataset ----------------------------------------------------------

win_res <- (win_import[,11:20])

win_res_trans <- t(apply(win_res, 1, function(x) {
  val = max(x)
  x[x != val] <- 0
  x[x == val] <- 1
  return(x)}
))

win_data <- data.frame(win_res_trans)
comb_data <- cbind(win_data, win_import[,21:31])

diff_win <- comb_data[comb_data$diff1_res == 1 | comb_data$diff2_res == 1 | comb_data$diff3_res == 1 | comb_data$diff4_res == 1 | 
            comb_data$diff5_res == 1 | comb_data$diff6_res == 1,]

parm <- list()
for(i in 1:10) {
  act_dist <- data.frame("parm" = unlist(diff_win[,11:21][i], use.names = FALSE),
                         "id" = paste0("Diff_Wins"),
                         "parm_name" = colnames(diff_win[11:21])[i])
  run_dist <- data.frame("parm" = runif(2000, low_parm[i], high_parm[i]),
                         "id" = paste0("Explored Space"),
                         "parm_name" = colnames(diff_win[11:21])[i])
  comb_dist <- rbind(act_dist, run_dist)
  parm[[i]] <- comb_dist
}



p_list <- list()

for(i in 1:10) {
  data_test <- parm[[i]]
  p_list[[i]] <- ggplot(data_test, aes(x=parm, fill=id)) + geom_density(alpha=.5) + theme_bw()  +
    theme(legend.text=element_text(size=10), axis.text.x=element_text(size=10),axis.ticks.y=element_blank(), axis.text.y=element_blank(),
          axis.title.y=element_text(size=10), axis.title.x= element_text(size=10), plot.margin = unit(c(0.25,0.4,0.15,0.55), "cm"),
          plot.title = element_text(size = 12, vjust = 3, hjust = 0.5, face = "bold")) + 
    scale_x_continuous(expand = c(0, 0), name = unique(data_test$parm_name)) + 
    scale_y_continuous(name = " ")
}

ggarrange(plotlist = p_list, nrow = 3, ncol = 4, common.legend = T, legend = "bottom")

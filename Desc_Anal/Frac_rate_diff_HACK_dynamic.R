library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis"); library("cowplot")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:6]),
                      sum(data_temp[4:6]),
                      sum(data_temp[3]),
                      sum(rowMeans(data_temp[4:6])),
                      sum(data_temp[4]),
                      sum(data_temp[5]),
                      sum(data_temp[6])),5)
  return(out_vec)
}

# ODE Equations -----------------------------------------------------------

#h_res <- eval(as.name(c("R1", "R2", "R3")[which.max(c(R1, R2, R3))]))
#m_res <- eval(as.name(c("R1", "R2", "R3")[setdiff(1:3, c(which.min(c(R1, R2, R3)), which.max(c(R1, R2, R3))))]))
#l_res <- eval(as.name(c("R1", "R2", "R3")[which.min(c(R1, R2, R3))]))

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    sigma_base1 <- sigma_1
    sigma_base2 <- sigma_2
    sigma_base3 <- sigma_3

    if(t > t_n) {
      sigma_1 <- sigma_base1*eff_tax1_1 
      sigma_2 <- sigma_base2*eff_tax2_1
      sigma_3 <- sigma_base3*eff_tax3_1
    }
    
    if(t > (t_n + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_2 
      sigma_2 <- sigma_base2*eff_tax2_2
      sigma_3 <- sigma_base3*eff_tax3_2
    }
    
    if(t > (t_n + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_3 
      sigma_2 <- sigma_base2*eff_tax2_3
      sigma_3 <- sigma_base3*eff_tax3_3
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_4 
      sigma_2 <- sigma_base2*eff_tax2_4
      sigma_3 <- sigma_base3*eff_tax3_4
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_5 
      sigma_2 <- sigma_base2*eff_tax2_5
      sigma_3 <- sigma_base3*eff_tax3_5
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      sigma_1 <- sigma_base1*eff_tax1_6 
      sigma_2 <- sigma_base2*eff_tax2_6
      sigma_3 <- sigma_base3*eff_tax3_6
    }
    
    dX = - beta*X*(W+(R1*c1)+(R2*c2)+(R3*c3)) + (1-sigma_1+sigma_2+sigma_3)*mu_w*W + (sigma_1+sigma_2+sigma_3)*mu_t*W +
      (1-sigma_2+sigma_3)*mu_r*R1 + (sigma_2+sigma_3)*mu_t*R1 +
      (1-sigma_1+sigma_3)*mu_r*R2 + (sigma_1+sigma_3)*mu_t*R2 + 
      (1-sigma_1+sigma_2)*mu_r*R3 + (sigma_1+sigma_2)*mu_t*R3
    
    dW = beta*X*W - (1-sigma_1+sigma_2+sigma_3)*mu_w*W - (sigma_1+sigma_2+sigma_3)*mu_t*W + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c1*R1 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c2*R2 + 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c3*R3 - 
      sigma_1*eta_wr*W - sigma_2*eta_wr*W - sigma_3*eta_wr*W
    
    dR1 = beta*X*R1*c1 - (1-sigma_2+sigma_3)*mu_r*R1 - (sigma_2+sigma_3)*mu_t*R1 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c1*R1 + sigma_1*eta_wr*W
    
    dR2 = beta*X*R2*c2 - (1-sigma_1+sigma_3)*mu_r*R2 - (sigma_1+sigma_3)*mu_t*R2 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c2*R2 + sigma_2*eta_wr*W
    
    dR3 = beta*X*R3*c3 - (1-sigma_1+sigma_2)*mu_r*R3 - (sigma_1+sigma_2)*mu_t*R3 - 
      (1-sigma_1+sigma_2+sigma_3)*eta_rw*c3*R3 + sigma_3*eta_wr*W
    
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX,dW,dR1,dR2,dR3)))
  })
}


# Dual Model --------------------------------------------------------------
init <- c(X = 0.99, W = 1-0.99, R1 = 0, R2 = 0, R3 = 0)

parms = c(beta = 5, sigma_1 = 0.25, sigma_2 = 0.25,sigma_3 = 0.25,
          mu_w = 1/12, mu_r = 1/10, 
          mu_t = 1/7, eta_wr = 0.01, eta_rw= 0.01, 
          c1= 0.95, c2 = 0.94, c3 = 0.9,
          eff_tax1_1 = 0.75, eff_tax2_1 = 0.75, eff_tax3_1 = 0.75, 
          eff_tax1_2 = 0.75, eff_tax2_2 = 0.75, eff_tax3_2 = 0.75, 
          eff_tax1_3 = 0.75, eff_tax2_3 = 0.75, eff_tax3_3 = 0.75, 
          eff_tax1_4 = 0.75, eff_tax2_4 = 0.75, eff_tax3_4 = 0.75, 
          eff_tax1_5 = 0.75, eff_tax2_5 = 0.75, eff_tax3_5 = 0.75, 
          eff_tax1_6 = 0.75, eff_tax2_6 = 0.75, eff_tax3_6 = 0.75, 
          PED1 = 1.5, 
          t_n = 3000)

multi_int_fun <- function(int_gen, basetax, time_between, parms, init, func){

  parms["time_between"] <- time_between
    
    #First Run
    testrun <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]]), parms = parms))
    values <- tail(testrun, 1)[4:6]
    
    print(values)
    
    parms[grep("eff_tax1", names(parms), value = T)] <- 1-(basetax*(values[1]/values[2]))
    parms[grep("eff_tax2", names(parms), value = T)] <- 1-(basetax*(values[2]/values[2]))
    parms[grep("eff_tax3", names(parms), value = T)] <- 1-(basetax*(values[3]/values[2]))
    parms[parms < 0] <- 0

    testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    
    #Second Run
    if(int_gen >= 2){
      testrun1 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + parms[["time_between"]]), parms = parms))
      values1 <- tail(testrun1, 1)[4:6]
      
      low_char <- names(values1)[which.min(values1)]
      high_char <- names(values1)[which.max(values1)]
      med_char <- names(values1)[setdiff(1:3, c(which.min(values1), which.max(values1)))]

      parms[grep(paste0("eff_tax",substr(low_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[low_char],1)/tail(testrun1[med_char],1)))
      parms[grep(paste0("eff_tax",substr(med_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[med_char],1)/tail(testrun1[med_char],1)))
      parms[grep(paste0("eff_tax",substr(high_char, 2, 2)), names(parms), value = T)[-1]] <- 1-(basetax*(tail(testrun1[high_char],1)/tail(testrun1[med_char],1)))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 3){
      testrun2 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*2)), parms = parms))
      values2 <- tail(testrun2, 1)[4:6]
      
      low_char1 <- names(values2)[which.min(values2)]
      high_char1 <- names(values2)[which.max(values2)]
      med_char1 <- names(values2)[setdiff(1:3, c(which.min(values2), which.max(values2)))]
      
      parms[grep(paste0("eff_tax",substr(low_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[low_char1],1)/tail(testrun2[med_char1],1)))
      parms[grep(paste0("eff_tax",substr(med_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[med_char1],1)/tail(testrun2[med_char1],1)))
      parms[grep(paste0("eff_tax",substr(high_char1, 2, 2)), names(parms), value = T)[-c(1,2)]] <- 1-(basetax*(tail(testrun2[high_char1],1)/tail(testrun2[med_char1],1)))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 4){
      testrun3 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*3)), parms = parms))
      values3 <- tail(testrun3, 1)[4:6]
      
      low_char2 <- names(values3)[which.min(values3)]
      high_char2 <- names(values3)[which.max(values3)]
      med_char2 <- names(values3)[setdiff(1:3, c(which.min(values3), which.max(values3)))]

      parms[grep(paste0("eff_tax",substr(low_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[low_char2],1)/tail(testrun3[med_char2],1)))
      parms[grep(paste0("eff_tax",substr(med_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[med_char2],1)/tail(testrun3[med_char2],1)))
      parms[grep(paste0("eff_tax",substr(high_char2, 2, 2)), names(parms), value = T)[-c(1,2,3)]] <- 1-(basetax*(tail(testrun3[high_char2],1)/tail(testrun3[med_char2],1)))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 5){
      testrun4 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*4)), parms = parms))
      values4 <- tail(testrun4, 1)[4:6]
      
      low_char3 <- names(values4)[which.min(values4)]
      high_char3 <- names(values4)[which.max(values4)]
      med_char3 <- names(values4)[setdiff(1:3, c(which.min(values4), which.max(values4)))]
      
      parms[grep(paste0("eff_tax",substr(low_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[low_char3],1)/tail(testrun4[med_char3],1)))
      parms[grep(paste0("eff_tax",substr(med_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[med_char3],1)/tail(testrun4[med_char3],1)))
      parms[grep(paste0("eff_tax",substr(high_char3, 2, 2)), names(parms), value = T)[-c(1,2,3,4)]] <- 1-(basetax*(tail(testrun4[high_char3],1)/tail(testrun4[med_char3],1)))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    if(int_gen >= 6){
      testrun5 <- data.frame(ode(y = init, func = func, times = seq(0, parms[["t_n"]] + (parms[["time_between"]]*5)), parms = parms))
      values5 <- tail(testrun5, 1)[4:6]
      
      low_char4 <- names(values5)[which.min(values5)]
      high_char4 <- names(values5)[which.max(values5)]
      med_char4 <- names(values5)[setdiff(1:3, c(which.min(values5), which.max(values5)))]
      
      parms[grep(paste0("eff_tax",substr(low_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[low_char4],1)/tail(testrun5[med_char4],1)))
      parms[grep(paste0("eff_tax",substr(med_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[med_char4],1)/tail(testrun5[med_char4],1)))
      parms[grep(paste0("eff_tax",substr(high_char4, 2, 2)), names(parms), value = T)[-c(1,2,3,4,5)]] <- 1-(basetax*(tail(testrun5[high_char4],1)/tail(testrun5[med_char4],1)))
      parms[parms < 0] <- 0
      
      testrun <- data.frame(ode(y = init, func = func, times = seq(0, 10000), parms = parms))
    }
    
    
    return(list(testrun, parms))
}


# Run the Differential Tax Model ------------------------------------------

test1 <- multi_int_fun(6, 0.5, 365*3, parms, init, amr)

#Calculate the Average

test <- test1[[1]]
write.csv(test, file = "/Users/amorgan/Documents/PostDoc/PrelimAnalysis/Theoretical_Analysis/Model_Output/dynamicoutput.csv")

test$AvgRes <- rowMeans(test[4:6])
test$TotInf <- rowSums(test[3:6])

comb_prop = melt(test, id.vars = "time", measure.vars = c("W","R1","R2","R3"))
trajdiff <- ggplot(data = comb_prop, aes(x = time, y = value*100, color = variable)) + geom_line() + theme_bw() +
  annotate("rect", xmin=c(parms[["t_n"]], parms[["t_n"]] + (365*3), parms[["t_n"]] + (365*3)*2, parms[["t_n"]] + (365*3)*3, parms[["t_n"]] + (365*3)*4, parms[["t_n"]] + (365*3)*5), 
           xmax= c(parms[["t_n"]] + (365*3), parms[["t_n"]] + (365*3)*2, parms[["t_n"]] + (365*3)*3, parms[["t_n"]] + (365*3)*4, parms[["t_n"]] + (365*3)*5, Inf), 
           ymin=c(rep(0, 6)), ymax=c(rep(Inf, 6)), alpha=rep(c(0.2, 0.125), 3), fill=rep("red", 6)) + 
  scale_x_continuous(name = "Time (days)", expand = c(0, 0), limits = c(0, 10000)) +  scale_y_continuous(name = "% Infected", expand = c(0, 0), limits = c(0,100)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + 
  scale_color_manual(values = c("black", "red", "blue", "green"))

diff_parms <- test1[[2]]

tax_data <- data.frame("time" = seq(0, 10000),
                       "R1" = c(rep(1, 3001), c(rep(diff_parms[["eff_tax1_1"]], 365*3), 
                                                        rep(diff_parms[["eff_tax1_2"]], 365*3),
                                                        rep(diff_parms[["eff_tax1_3"]], 365*3),
                                                        rep(diff_parms[["eff_tax1_4"]], 365*3),
                                                        rep(diff_parms[["eff_tax1_5"]], 365*3),
                                                        rep(diff_parms[["eff_tax1_6"]], (365*3)+430))), 
                       "R2" = c(rep(1, 3001), c(rep(diff_parms[["eff_tax2_1"]], 365*3), 
                                                        rep(diff_parms[["eff_tax2_2"]], 365*3),
                                                        rep(diff_parms[["eff_tax2_3"]], 365*3),
                                                        rep(diff_parms[["eff_tax2_4"]], 365*3),
                                                        rep(diff_parms[["eff_tax2_5"]], 365*3),
                                                        rep(diff_parms[["eff_tax2_6"]], (365*3)+430))), 
                       "R3" = c(rep(1, 3001), c(rep(diff_parms[["eff_tax3_1"]], 365*3), 
                                                        rep(diff_parms[["eff_tax3_2"]], 365*3),
                                                        rep(diff_parms[["eff_tax3_3"]], 365*3),
                                                        rep(diff_parms[["eff_tax3_4"]], 365*3),
                                                        rep(diff_parms[["eff_tax3_5"]], 365*3),
                                                        rep(diff_parms[["eff_tax3_6"]], (365*3)+430))))

tax_melt <- melt(tax_data, id.vars = c("time"), measure.vars = c("R1", "R2", "R3"))
taxdiff_comb <- ggplot(tax_melt, aes(x = time, y = (1-value)*100, color = variable)) + geom_line(size = 1.2) + theme_bw() + 
  annotate("rect", xmin=c(parms[["t_n"]], parms[["t_n"]] + (365*3), parms[["t_n"]] + (365*3)*2, parms[["t_n"]] + (365*3)*3, parms[["t_n"]] + (365*3)*4, parms[["t_n"]] + (365*3)*5), 
           xmax= c(parms[["t_n"]] + (365*3), parms[["t_n"]] + (365*3)*2, parms[["t_n"]] + (365*3)*3, parms[["t_n"]] + (365*3)*4, parms[["t_n"]] + (365*3)*5, Inf), 
           ymin=c(rep(0, 6)), ymax=c(rep(Inf, 6)), alpha=rep(c(0.2, 0.125), 3), fill=rep("red", 6)) + 
  scale_x_continuous(name = "Time (days)", expand = c(0, 0)) +  scale_y_continuous(name = "Effective Tax Rate (%)", expand = c(0, 0), limits = c(0,101)) + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x = element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "bottom") + 
  scale_color_manual(values = c("red", "blue", "green"))

ggarrange(trajdiff, taxdiff_comb, nrow = 2, ncol = 1, common.legend = T, legend = "bottom")

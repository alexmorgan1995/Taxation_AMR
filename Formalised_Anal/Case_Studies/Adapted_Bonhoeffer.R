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

# 3 Strain Bonhoeffer ODEs ------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    f1_base <- f_1
    f2_base <- f_2
    f3_base <- f_3
    
    if(t > t_n) {
      f_1 <- f1_base*(1-(eff_tax1_1*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_1*PED2))
      f_3 <- f3_base*(1-(eff_tax3_1*PED3))
    }
    
    if(t > (t_n + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_2*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_2*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_2*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_3*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_3*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_3*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_4*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_4*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_4*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_5*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_5*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_5*PED3)) 
    }
    
    if(t > (t_n + time_between + time_between + time_between + time_between + time_between)) {
      f_1 <- f1_base*(1-(eff_tax1_6*PED1)) 
      f_2 <- f2_base*(1-(eff_tax2_6*PED2)) 
      f_3 <- f3_base*(1-(eff_tax3_6*PED3)) 
    }

    f_1 <- ifelse(f_1 > 0, f_1, 0)
    f_2 <- ifelse(f_2 > 0, f_2, 0)
    f_3 <- ifelse(f_3 > 0, f_3, 0)
    
    dx = lambda - lambda*x - beta*x*(y_w + y_1 + y_2 + y_3 + y_12 + y_13 + y_23 + y_123) + 
      r_w*y_w + r_1*y_1 + r_2*y_2 + r_3*y_3 + r_12*y_12 + r_23*y_23 + r_13*y_13 + r_123*y_123 +
      h*(1-s)*((f_1+f_2+f_3)*y_w + f_1*(y_2+y_3+y_23)+ f_2*(y_1+y_3+y_13) + f_3*(y_1+y_2+y_12))
    
    dy_w = (beta*x - lambda - r_w - h*(f_1 + f_2 + f_3))*y_w
    
    dy_1 = (beta*x - lambda - r_1 - h*(f_2 + f_3))*y_1 + h*s*f_1*y_w
    dy_2 = (beta*x - lambda - r_2 - h*(f_1 + f_3))*y_2 + h*s*f_2*y_w
    dy_3 = (beta*x - lambda - r_3 - h*(f_2 + f_3))*y_3 + h*s*f_3*y_w
    
    dy_12 = (beta*x - lambda - r_12 - h*f_3)*y_12 + h*s*(f_1*y_2 + f_2*y_1)
    dy_13 = (beta*x - lambda - r_13 - h*f_2)*y_13 + h*s*(f_1*y_3 + f_3*y_1)
    dy_23 = (beta*x - lambda - r_23 - h*f_1)*y_23 + h*s*(f_3*y_2 + f_2*y_3)
    
    dy_123 = (beta*x - lambda - r_123)*y_123 + h*s*(f_1*y_23 + f_2*y_13 + f_3*y_12)
    
    #Calculating the Proportion Integrals
    
    return(list(c(dx,dy_w,dy_1,dy_2,dy_3, dy_12, dy_13, dy_23, dy_123)))
  })
}

# Baseline Parms ----------------------------------------------------------

init <- c(x = 0.99, y_w = 1-0.99, y_1 = 0, y_2 = 0, y_3 = 0, y_12 = 0, y_13 = 0, y_23 = 0, y_123 = 0)

parms = c(beta = 1, f_1 = 0.25, f_2 = 0.25, f_3 = 0.25,
          lambda = 0.5, h = 0.1, s = 0.8,
          r_w = 1/13, 
          r_1 = 1/12, r_2 = 1/11.5, r_3 = 1/11,
          r_12 = 1/10, r_13 = 1/9.5, r_23 = 1/9,
          r_123 = 1/2,
          eff_tax1_1 = 0, eff_tax2_1 = 0, eff_tax3_1 = 0, 
          eff_tax1_2 = 0, eff_tax2_2 = 0, eff_tax3_2 = 0, 
          eff_tax1_3 = 0, eff_tax2_3 = 0, eff_tax3_3 = 0, 
          eff_tax1_4 = 0, eff_tax2_4 = 0, eff_tax3_4 = 0, 
          eff_tax1_5 = 0, eff_tax2_5 = 0, eff_tax3_5 = 0, 
          eff_tax1_6 = 0, eff_tax2_6 = 0, eff_tax3_6 = 0, 
          PED1 = 1, PED2 = 1, PED3 = 1, 
          t_n = 3000, time_between = Inf)

# Run the Models ----------------------------------------------------------

#Flat Tax
parms[grep("eff_tax", names(parms), value = TRUE)] <- 0.5
testrun_flat <- data.frame(ode(y = init, func = amr, times = seq(0, 10000), parms = parms))

melt_t <- melt(testrun_flat, measure.vars = colnames(testrun_flat)[-1], id.vars = "time")
#melt_t <- melt(testrun_flat, measure.vars = c("y_w", "y_1", "y_2", "y_3"), id.vars = "time")

ggplot(melt_t, aes(time, value, color = variable)) + theme_bw() + geom_line()

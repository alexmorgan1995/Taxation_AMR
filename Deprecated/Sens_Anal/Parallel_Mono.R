library("deSolve"); library("ggplot2"); library("reshape2"); library("ggpubr"); library("rootSolve"); library("viridis")
library("sensitivity"); library("parallel")

rm(list=ls())

setwd("/Users/amorgan/Documents/PostDoc/PrelimAnalysis/RCode")

# ODE Equations -----------------------------------------------------------

amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    if(t > t_n) {
      sigma_1 <- sigma_1*eff_tax1  
      sigma_2 <- sigma_2*eff_tax2
    }
    
    dX = - beta*X*(W+(R1*c1)+(R2*c2)) + (1-sigma_1+sigma_2)*mu_w*W + (sigma_1+sigma_2)*mu_t*W +
      (1-sigma_2)*mu_r*R1 + sigma_2*mu_t*R1 +
      (1-sigma_1)*mu_r*R2 + sigma_1*mu_t*R2
    
    dW = beta*X*W - (1-sigma_1+sigma_2)*mu_w*W - (sigma_1+sigma_2)*mu_t*W + 
      (1-sigma_1+sigma_2)*eta_rw*c1*R1 + (1-sigma_1+sigma_2)*eta_rw*c2*R2 - 
      sigma_1*eta_wr*W - sigma_2*eta_wr*W
    
    dR1 = beta*X*R1*c1 - (1-sigma_2)*mu_r*R1 - sigma_2*mu_t*R1 - (1-sigma_1+sigma_2)*eta_rw*c1*R1 + sigma_1*eta_wr*W
    
    dR2 = beta*X*R2*c2 - (1-sigma_1)*mu_r*R2 - sigma_1*mu_t*R2 - (1-sigma_1+sigma_2)*eta_rw*c2*R2 + sigma_2*eta_wr*W
    
    #Calculating the Proportion Integrals
    
    return(list(c(dX,dW,dR1,dR2)))
  })
}

# Integral Function -------------------------------------------------------

integral <- function(data, t_n){
  data_temp <- data[data[,1] > t_n,]
  
  out_vec <- signif(c(sum(data_temp[3:5]),
                      sum(rowMeans(data_temp[4:5])),
                      sum(data_temp[4]),
                      sum(data_temp[5])),5)
  return(out_vec)
}

# Parameter Set -----------------------------------------------------------

parm_vector_func <- function(n, parms, min_max) {
  
  reord_parm <- data.frame(matrix(nrow = 0, ncol = length(parms)))
  
  for(i in 1:length(min_max)) {
    parm_name <- names(parms)[i]
    oth_parm <- data.frame(t(matrix(parms[names(parms) != parm_name])))[rep(1,n+1), ]
    parm_rep <- data.frame(parm_name = seq(from = min_max[[i]][1], to = min_max[[i]][2], 
                                     by = (min_max[[i]][2]-min_max[[i]][1])/n))
    comb <- cbind(oth_parm, parm_rep)
    
    colnames(comb) <- c(names(parms)[names(parms) != parm_name], parm_name)
    
    temp_out <- comb[,names(parms)]
    reord_parm <- rbind(reord_parm, temp_out)
  } 
  return(reord_parm)
}

out <- parm_vector_func(1000, 
                        
                        parms = c(beta = 1, sigma_1 = 0.35, sigma_2 = 0.35, mu_w = 1/15, mu_r = 1/10, 
                                  mu_t = 1/5, eta_wr = 0.02, eta_rw= 0.02, c1= 0.9, c2 = 0.85,
                                  eff_tax1 = 0.5, eff_tax2 = 0.5, t_n = 1000),
                        
                        min_max = list(c(0,10),
                                       c(0,1),
                                       c(0,1),
                                       c(1/50,1/0.5),
                                       c(1/50,1/0.5),
                                       c(1/50,1/0.5),
                                       c(0.002,0.2),
                                       c(0.002,0.2),
                                       c(0,1),
                                       c(0,1),
                                       c(0,1),
                                       c(0,1)))

# Creating the Parallel Montonicity Function ------------------------------

mono_func <- function(n, parms_frame, init) {
    parms = parms_frame[n,]
    run <- data.frame(ode(y = init, func = amr, times = seq(0, 2000), parms = parms))
    return(c(integral(run, parms[["t_n"]]), parms_frame[n,]
             ))
}

# Running the Parallelisable Function -------------------------------------

output <- mclapply(1:nrow(out), mono_func, parms_frame = out, init = c(X = 0.99, W = 1-0.99, R1 = 0, R2 =0),
                   mc.cores = 6)

fin_data <- as.data.frame(cbind(do.call("rbind", output), as.character(rep(c("beta", "sigma_1", "sigma_2", "mu_w", "mu_r", 
                                                                "mu_t", "eta_wr", "eta_rw", "c1", "c2",
                                                                "eff_tax1", "eff_tax2"), rep(1001, 12)))))

fin_data <- as.data.frame(lapply(fin_data,unlist))
colnames(fin_data)[1:4] <- c("totinf", "avg", "R1", "R2")
colnames(fin_data)[18] <- "parm"

na <- fin_data[rowSums(is.na(fin_data))==0,]

# Plotting ----------------------------------------------------------------

p_list <- list()

for(i in 1:4) {
  sub_list <- list()
  data <- cbind(fin_data[i], fin_data[5:18])
  
  for(j in 1:length(unique(fin_data$parm))){
  parm <- unique(fin_data$parm)[j]
   test <- data.frame("outcome" = data[data$parm == parm,][,1],
               "parm" = data[data$parm == parm,][[parm]])
   
   sub_list[[j]] <- ggplot(test, aes(x = parm, y = outcome)) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
     scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, max(na.omit(test$outcome)*1.1)), expand = c(0, 0)) +
     labs(x = parm, y = names(fin_data)[i]) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"))
  }
  p_list[[i]] <- sub_list
}

# Plotting Combined -------------------------------------------------------

gg_list <- list()

for(i in 1:4){
  gg_list[[i]] <- ggarrange(plotlist = p_list[[i]], nrow = 4, ncol = 3)
}

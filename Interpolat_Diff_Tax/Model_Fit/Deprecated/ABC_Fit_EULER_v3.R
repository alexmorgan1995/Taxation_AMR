library("deSolve"); library("bayestestR"); library("tmvtnorm"); library("parallel")
rm(list=ls())

# Model -------------------------------------------------------------------

amr <- function(t, y, parms, sigma_use1, sigma_use2, sigma_use3) {
  with(as.list(c(y, parms)), {
    
    #Specify the time-varying functions
    
    sigma_use1 <- sigma_func1(t)
    sigma_use2 <- sigma_func2(t)
    sigma_use3 <- sigma_func3(t)
    
    #ODES Below
    
    dX = lambda - lambda*X - beta*X*(Wt + R1*c1 + R2*c2 + R3*c3 + R12*c12 + R13*c13 + R23*c23 + R123*c123) +
      r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + r_r*R1*(1-(sigma_use2 + sigma_use3)) + r_r*R2*(1-(sigma_use1 + sigma_use3)) + r_r*R3*(1-(sigma_use1 + sigma_use2)) +
      r_rr*R12*(1-sigma_use3) + r_rr*R13*(1-sigma_use2) + r_rr*R23*(1-sigma_use1) + 
      r_rrr*R123 + 
      r_t*(1-rho)*(Wt*(sigma_use1 + sigma_use2 + sigma_use3) + R1*(sigma_use2 + sigma_use3) + 
                     R2*(sigma_use1 + sigma_use3) + R3*(sigma_use1 + sigma_use2) + 
                     R12*sigma_use3 + R13*sigma_use2 + R23*sigma_use1)
    
    dWt = - lambda*Wt + beta*X*Wt - r_wt*Wt*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) - r_t*Wt*(1-rho)*(sigma_use1 + sigma_use2 + sigma_use3) +
      eta_rw*(R1 + R2 + R3 + R12 + R13 + R23 + R123)*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) - 
      eta_wr*Wt*rho*(sigma_use1 + sigma_use2 + sigma_use3)
    
    dR1 = - lambda*R1 + beta*X*R1*c1 - r_t*(1-rho)*(sigma_use2 + sigma_use3)*R1 - r_r*(1-(sigma_use2 + sigma_use3))*R1 - eta_rr*R1*rho*sigma_use2 -
      eta_rr*R1*rho*sigma_use3 - eta_rw*R1*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use1
    
    dR2 = - lambda*R2 + beta*X*R2*c2 - r_t*(1-rho)*(sigma_use1 + sigma_use3)*R2 - r_r*(1-(sigma_use1 + sigma_use3))*R2 - eta_rr*R2*rho*sigma_use1 -
      eta_rr*R2*rho*sigma_use3 - eta_rw*R2*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use2
    
    dR3 = - lambda*R3 + beta*X*R3*c3 - r_t*(1-rho)*(sigma_use1 + sigma_use2)*R3 - r_r*(1-(sigma_use1 + sigma_use2))*R3 - eta_rr*R3*rho*sigma_use1 -
      eta_rr*R3*rho*sigma_use2 - eta_rw*R3*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_wr*rho*Wt*sigma_use3
    
    
    dR12 = - lambda*R12 + beta*X*R12*c12 - r_t*(1-rho)*sigma_use3*R12 - r_rr*(1-sigma_use3)*R12 - eta_rrr*R12*rho*sigma_use3 -
      eta_rw*R12*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R1*rho*sigma_use2 + eta_rr*R2*rho*sigma_use1 
    
    dR13 = - lambda*R13 + beta*X*R13*c13 - r_t*(1-rho)*sigma_use2*R13 - r_rr*(1-sigma_use2)*R13 - eta_rrr*R13*rho*sigma_use2 -
      eta_rw*R13*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R1*rho*sigma_use3 + eta_rr*R3*rho*sigma_use1 
    
    dR23 = - lambda*R23 + beta*X*R23*c23 - r_t*(1-rho)*sigma_use1*R23 - r_rr*(1-sigma_use1)*R23 - eta_rrr*R23*rho*sigma_use1 -
      eta_rw*R23*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + eta_rr*R2*rho*sigma_use3 + eta_rr*R3*rho*sigma_use2 
    
    dR123 = - lambda*R123 + beta*X*R123*c123 - r_rrr*R123 - eta_rw*R123*(1 - (sigma_use1 + sigma_use2 + sigma_use3)) + 
      eta_rrr*rho*(sigma_use3*R12 + sigma_use2*R13 + sigma_use1*R23)
    
    return(list(c(dX,dWt,
                  dR1,dR2,dR3,
                  dR12,dR13,dR23,
                  dR123)))
  })
}

# Function to Aggregate Resistance ----------------------------------------

agg_func <- function(data) {
  agg_data <- data.frame("time" = data$time,
                         "Susc" = data$X,
                         "WT" = data$Wt, 
                         "R1" = data$R1 + data$R12 + data$R13 + data$R123,
                         "R2" = data$R2 + data$R12 + data$R23 + data$R123,
                         "R3" = data$R3 + data$R13 + data$R23 + data$R123)
  return(agg_data)
}

# Extract Sigmas for the ApproxFun Function -------------------------------

approx_sigma <- function(sigma_mat){
  
  usage = data.frame("time" = seq(0,10000),
                     "PopUsage1" = c(rep(sigma_mat[1,1], 3000),
                                     rep(sigma_mat[1,2], 365*3), rep(sigma_mat[1,3], 365*3), rep(sigma_mat[1,4], 365*3),
                                     rep(sigma_mat[1,5], 365*3), rep(sigma_mat[1,6], 365*3), rep(sigma_mat[1,7], 10001 - (3000 + (365*3)*5))),
                     
                     "PopUsage2" = c(rep(sigma_mat[2,1], 3000),
                                     rep(sigma_mat[2,2], 365*3), rep(sigma_mat[2,3], 365*3), rep(sigma_mat[2,4], 365*3),
                                     rep(sigma_mat[2,5], 365*3), rep(sigma_mat[2,6], 365*3), rep(sigma_mat[2,7], 10001 - (3000 + (365*3)*5))),
                     
                     "PopUsage3" = c(rep(sigma_mat[3,1], 3000),
                                     rep(sigma_mat[3,2], 365*3), rep(sigma_mat[3,3], 365*3), rep(sigma_mat[3,4], 365*3),
                                     rep(sigma_mat[3,5], 365*3), rep(sigma_mat[3,6], 365*3), rep(sigma_mat[3,7], 10001 - (3000 + (365*3)*5))))
  return(usage)
}

# ODE Wrapper Function ----------------------------------------------------

ode_wrapper <- function(times, y, parms, func, approx_sigma) {
  
  sigma_mat = matrix(c(rep(parms[["sigma1"]], 7),
                       rep(parms[["sigma2"]], 7),
                       rep(parms[["sigma3"]], 7)), 
                     nrow = 3, ncol = 7, byrow = T)
  eff_tax <- parms[["eff_tax"]]; PED <- parms[["PED"]]
  
  if(parms[["int_round"]] > 0 ) {
    for(i in 1:parms[["int_round"]]) {
      stor_sigma <- sigma_mat[,i]
      
      sigma_mat[,(i+1):7] = c(stor_sigma[1]*(1 + ((eff_tax[1,i]*PED[1,1]) + (eff_tax[2,i]*PED[2,1]) + (eff_tax[3,i]*PED[3,1]))),
                              stor_sigma[2]*(1 + ((eff_tax[1,i]*PED[1,2]) + (eff_tax[2,i]*PED[2,2]) + (eff_tax[3,i]*PED[3,2]))),
                              stor_sigma[3]*(1 + ((eff_tax[1,i]*PED[1,3]) + (eff_tax[2,i]*PED[2,3]) + (eff_tax[3,i]*PED[3,3]))))

      sigma_mat[,(i+1):7][sigma_mat[,(i+1)] < 0.01] <- 0.01
      
      if(colSums(sigma_mat)[i+1] > 1) {
        sigma_mat[,(i+1):7] <- sigma_mat[,i+1]/(sum(sigma_mat[,i+1])+0.01)
      }
    }
  }
  
  parms[["sigma_mat"]] <- sigma_mat
  
  sigma_data <- approx_sigma(sigma_mat)
  
  sigma_func1 <<- approxfun(sigma_data[,c(1,2)], rule = 2)
  sigma_func2 <<- approxfun(sigma_data[,c(1,3)], rule = 2)
  sigma_func3 <<- approxfun(sigma_data[,c(1,4)], rule = 2)
  
  #Run the model 
  out <- data.frame(ode(y = y, func = func, times = times, parms = parms))
  
  return(list(out, parms))
}

# Function to Aggregate Resistance ----------------------------------------

#Compute the distances for all 3 summary statistics - this section involves running the model

computeDistanceABC_ALEX <- function(fitmodel, ode_wrapper, approx_sigma, thetaparm, init.state) {
  out <- ode_wrapper(y = init.state, func = fitmodel, times = seq(0, 10000), parms = thetaparm, approx_sigma)[[1]]
  end_res <- as.numeric(tail(agg_func(out)[4:6], 1))
  return(abs(c(end_res[1] - 0.4,
           end_res[2] - 0.25,
           end_res[3] - 0.1)))
}

#Where G is the number of generations
#Function to 100% make sure the sampled particles for all parameters are non zero
prior.non.zero<-function(par, lm.low, lm.upp){
  prod(sapply(1:11, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

# Single Particle of the Model Fit ---------------------------------------------

singlerun <- function(x, G, init.state, distanceABC, fitmodel, thetaparm, epsilon, 
                      lm.low, lm.upp, w.old, sigma, res.old, N, ode_wrapper, approx_sigma) {
  i <- 0
  m <- 0
  w.new <- 0
  
  while(i <= 1) {
    
    m <- m + 1
    
    if(G == 1) {
      etas <- runif(2, min = 0, max = 4)
      cs <- rev(sort(runif(7, 0.5, 1)))
      
      d_beta <- runif(1, min = 0, max = 10)
      d_eta_wr <- max(etas)
      d_eta_rw <- min(etas)
      d_eta_rr_rrr <- runif(1, 0, 0.5)
      d_c1 <- cs[1]
      d_c2 <- cs[2]
      d_c3 <- cs[3]
      d_c12 <- cs[4]
      d_c13 <- cs[5]
      d_c23 <- cs[6]
      d_c123 <- cs[7]
      
    } else { 
      p <- sample(seq(1,N),1,prob= w.old) # check w.old here
      par <- rtmvnorm(1, mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
      
      #at this point rearrange the parameters 
      sorted_pars <- rev(sort(c(par[5], par[6], par[7], par[8], par[9], par[10], par[11])))
      d_beta <- par[1]
      d_eta_wr <- max(c(par[2],par[3]))
      d_eta_rw <- min(c(par[2],par[3]))
      d_eta_rr_rrr <- par[4]
      d_c1 <- sorted_pars[1]
      d_c2 <- sorted_pars[2]
      d_c3 <- sorted_pars[3]
      d_c12 <- sorted_pars[4]
      d_c13 <- sorted_pars[5]
      d_c23 <- sorted_pars[6]
      d_c123 <- sorted_pars[7]
    }
    
    new.parms = c(d_beta, d_eta_wr, d_eta_rw, d_eta_rr_rrr, d_c1, d_c2, d_c3, d_c12, d_c13, d_c23, d_c123)
    
    if(prior.non.zero(new.parms, lm.low, lm.upp)) {
      
      thetaparm["beta"] <- new.parms[1]
      thetaparm["eta_wr"] <- new.parms[2]
      thetaparm["eta_rw"] <- new.parms[3]
      thetaparm["eta_rr"] <- new.parms[4]
      thetaparm["eta_rrr"] <- new.parms[4]
      thetaparm["c1"] <- new.parms[5]
      thetaparm["c2"] <- new.parms[6]
      thetaparm["c3"] <- new.parms[7]
      thetaparm["c12"] <- new.parms[8]
      thetaparm["c13"] <- new.parms[9]
      thetaparm["c23"] <- new.parms[10]
      thetaparm["c123"] <- new.parms[11]
      
      dist_mod <- distanceABC(fitmodel, ode_wrapper, approx_sigma, thetaparm, init.state)
      
      if((dist_mod[1] <= epsilon[["R1"]][G]) && (dist_mod[2] <= epsilon[["R2"]][G]) && (dist_mod[3] <= epsilon[["R3"]][G]) && (!is.na(dist_mod))) {
        if(G==1){
          w.new <- 1
        } else {
          w1 <- prod(sapply(c(1:11), function(b) dunif(new.parms[b], min=lm.low[b], max=lm.upp[b])))
          w2 <- sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(new.parms, mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
          w.new <- w1/w2
        }
        i <- i + 1
        return(list(dist_mod, m, new.parms, w.new))
      }
    }
  }
}

# ABC-SMC Function --------------------------------------------------------

ABC_algorithm <- function(N, G, single_run, distanceABC, fitmodel, init.state, epsilon, lm.low, lm.upp, thetaparm, ode_wrapper, approx_sigma)  {
  out <- list()
  
  for(g in 1:G) {
    
    print(paste0("Generation ", g, " | Time: ", Sys.time()))
    
    if(g == 1) {
      sigma <- 0
      res.old <- 0
      w.old <- 0
    }

    particles <- mclapply(1:N, 
                          G = g,
                          FUN = singlerun, 
                          distanceABC = distanceABC,
                          fitmodel = fitmodel, 
                          init.state = init.state, 
                          thetaparm = thetaparm, 
                          epsilon = epsilon,
                          lm.low = lm.low,
                          lm.upp = lm.upp,
                          w.old = w.old, 
                          sigma = sigma, 
                          res.old = res.old,
                          N = N, 
                          ode_wrapper = ode_wrapper, 
                          approx_sigma = approx_sigma,
                          mc.cores = 10)
    
    print(particles)
    dat_dist <- as.matrix(do.call(rbind, lapply(particles, "[[", 1)))
    dat_nruns <- do.call(sum, lapply(particles, "[[", 2))
    res.new <- as.matrix(do.call(rbind, lapply(particles, "[[", 3)))
    w.new <- as.matrix(do.call(rbind, lapply(particles, "[[", 4)))
    
    sigma <- cov(res.new) 
    res.old <- res.new
    w.old <- w.new/sum(w.new)
    
    out[[g]] <- list(dat_nruns, dat_dist, res.old, w.old)
    colnames(res.old) <- c("beta", "eta_wr", "eta_rw", "eta_rr_rrr", "c1", "c2", "c3", "c12", "c13", "c23", "c123")
    write.csv(res.old, file =  paste("/cluster/home/amorgan/Model_Fit_Output/ABC_v3_",g,".csv",sep=""), row.names=FALSE)
  }
  return(out)
}

# Running the Model Fit ---------------------------------------------------

ABC_algorithm(N = 1000,
              G = 10,
              single_run = single_run,
              distanceABC = computeDistanceABC_ALEX, 
              fitmodel = amr, 
              init.state = c(X = 0.99, Wt = 1-0.99, R1 = 0, R2 = 0, R3 = 0,
                             R12 = 0, R13 = 0, R23 = 0,
                             R123 = 0), 
              epsilon = list("R1" =  c(0.5*1, 0.5*0.8, 0.5*0.6, 0.5*0.4, 0.5*0.3, 0.5*0.2, 0.5*0.1, 0.5*0.05, 0.5*0.025, 0.5*0.01),
                             "R2" = c(0.25*1, 0.25*0.8, 0.25*0.6, 0.25*0.4, 0.25*0.3, 0.25*0.2, 0.25*0.1, 0.25*0.05, 0.25*0.025, 0.25*0.01),
                             "R3" =  c(0.1*1, 0.1*0.8, 0.1*0.6, 0.1*0.4, 0.1*0.3, 0.1*0.2, 0.1*0.1, 0.1*0.05, 0.1*0.025, 0.1*0.01)),
              lm.low = c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
              lm.upp = c(10, 4, 4, 0.5, 1, 1, 1, 1, 1, 1, 1),
              thetaparm = list(lambda = 1/365*(2), int_round = 1, 
                               sigma1 = 0.25, sigma2 = 0.25, sigma3 = 0.25,
                               r_wt = 1/12, r_r = 1/10,  r_rr = 1/9,  r_rrr = 1/8, 
                               r_t = 1/7,
                               PED = matrix(c(-1, 0.4, 0.4, 
                                              0.4, -1, 0.4,
                                              0.4, 0.4, -1), #Be aware of this matrix
                                            nrow = 3, ncol = 3, byrow = T),
                               eff_tax = matrix(c(0, 0, 0, 0, 0, 0, 
                                                  0, 0, 0, 0, 0, 0, 
                                                  0, 0, 0, 0, 0, 0), 
                                                nrow = 3, ncol = 6, byrow = T),
                               t_n = 3000, time_between = Inf, rho = 0.05, base_tax = 0.5), ode_wrapper, approx_sigma)

saveRDS(test, file = "/cluster/home/amorgan/Model_Fit_Output/dist_list_v3.rds")

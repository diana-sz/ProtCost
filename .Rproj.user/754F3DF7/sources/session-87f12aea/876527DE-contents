rm(list=ls(all=TRUE))

library(here)
# require('rstudioapi') 
# require('readODS')
# require('nloptr')
# require('lpSolve')

directory <- paste0(here(), "/code")
setwd(directory) 

is.reversible <- 1
predict.parameters <- 0
modelname <- "A8alt_trans"

source("initialize_model.R")
n_conditions <- 1
rho_cond <- rho_cond[1]

q0_wt <- q0
last_feasible_q0 <- q0


results_list <- list()

transporter <- which(reaction == "TC")
transporter2 <- which(reaction == "TC2")
kcat_tC1 <- kcatf[transporter]
km_tC1 <- KS[1,transporter]

phis_to_test <- c(0, 0.0001, seq(0.01, 0.5, 0.01))
kms_to_test <- c(km_tC1/10, km_tC1, km_tC1*10)
kcats_to_test <- c(kcat_tC1*10, kcat_tC1, kcat_tC1/10)

for (km in kms_to_test){
  KS[1,transporter2] <- km
  
  for (kcat_t in kcats_to_test){
    kcatf[transporter2] <- kcat_t
    kcatb[transporter2] <- kcat_t/5
    
    # get optimal phis
    source("solver_loop.R")
    
    p_opt <- prot(q_opt[1,])
    opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
    opt_phi <- opt_phis[transporter]

    for (fraction in phis_to_test){
      min_phi[transporter] <- fraction
      max_phi[transporter] <- fraction+1e-5
      
      source("solver_loop.R")
      
      fs <- q_opt[1, ]
      results_list[[length(results_list) + 1]] <- data.frame(
        kcat = kcat_t,
        km = km,
        phi = fraction,
        rel_phi = fraction/opt_phi,
        mu = mu_opt,
        convergence = res$convergence,
        t(c(
          setNames(fs, paste0("f.", reaction)),
          setNames(v(fs), paste0("v.", reaction)),
          setNames(prot(fs), paste0("p.", reaction)),
          setNames(c(a, ci(fs)), paste0("c.", reactant))
        ))
      )
      
    }
    
    # reset phis
    min_phi[] <- 0
    max_phi[] <- 1
  }
  
}

results <- do.call(rbind, results_list)
write.csv(results, paste0("../data/", modelname, "_alt_trans.csv"))


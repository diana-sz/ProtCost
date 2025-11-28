rm(list=ls(all=TRUE))

library(here)
# library(rstudioapi) 
# library(readODS)
# library(nloptr)

directory <- paste0(here(), "/code")
setwd(directory) 

phis_to_test <- c(seq(0, 0.8, 0.01), seq(0.8, 0, -0.01))
predict.parameters <- 0
primary_c_source <- 0.5

alt_concentrations <- c(0.001, 0.01, 0.05, 0.1, 0.2, 1, 10)


for(is.reversible in c(1,0)){
  results_list <- list()
  
  for (x_C2 in alt_concentrations){
    modelname <- "A8dekel"
    source("initialize_model.R")
    q0_wt <- q0
    last_feasible_q0 <- q0
    
    alt_ind <- which(reaction == "LAC")
    a_cond[1,1] <- primary_c_source
    a_cond[2,1] <- x_C2
    n_conditions <- 1
    
    source("GBA_solver.R")
    
    mu_orig <- mu_opt
    p_opt <- prot(q_opt[1,])
    opt_phis <- (p_opt/rho_cond)/sum(p_opt/rho_cond)
    opt_phi <- opt_phis[alt_ind]
    if(opt_phi < 1e-8){
      opt_phi_nonzero <- FALSE
      opt_phi <- 0.01
    }else{
      opt_phi_nonzero <- TRUE
    }
    
    for (fraction in phis_to_test){
      print(paste("LAC=", x_C2, "fraction=", fraction))
      
      min_phi[alt_ind] <- fraction
      max_phi[alt_ind] <- fraction+1e-5
      
      source("solver_loop.R")
      
      fs <- q_opt[1, ]
      results_list[[length(results_list) + 1]] <- data.frame(
        x_C = a_cond[1,1],
        x_C2 = x_C2,
        protein = "LAC",
        phi = fraction,
        rel_phi = fraction / opt_phi,
        mu = mu_opt,
        mu_norm = mu_opt / mu_orig,
        opt_phi_nonzero = opt_phi_nonzero,
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
  
  results <- do.call(rbind, results_list)
  
  write.csv(results, paste0("../data/", modelname, ".csv"))
  
}



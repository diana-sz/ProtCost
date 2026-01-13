rm(list=ls(all=TRUE))

library(here)

directory <- paste0(here(), "/code")
setwd(directory) 

predict.parameters <- 0
phis_to_test <- seq(0.1, 0.4, 0.1)
gfp_fracion <- seq(0, 1, 0.2)

for(is.reversible in c(1,0)){
  modelname <- "M10_GFP_RFP"
  
  source("initialize_model.R")
  n_conditions <- 1
  rho_cond <- rho_cond[1]
  
  q0_wt <- q0
  last_feasible_q0 <- q0
  
  results_list <- list()
  gfp_ind <- which(reaction == "GFP")
  rfp_ind <- which(reaction == "RFP")
  for(fraction in phis_to_test){
    
    for (gfp in gfp_fracion){
      min_phi[gfp_ind] <- gfp*fraction
      max_phi[gfp_ind] <- gfp*fraction+1e-5
      min_phi[rfp_ind] <- (1-gfp)*fraction
      max_phi[rfp_ind] <- (1-gfp)*fraction+1e-5
      
      source("solver_loop.R")
      
      fs <- q_opt[1, ]
      results_list[[length(results_list) + 1]] <- data.frame(
        gfp_phi = gfp*fraction,
        rfp_phi = (1-gfp)*fraction,
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
  }
  results <- do.call(rbind, results_list)
  write.csv(results, paste0("../data/", modelname, ".csv"))
}




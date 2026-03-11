# Test cost of an idle protein

rm(list=ls(all=TRUE))

library(here)

predict.parameters <- 0
phis_to_test <- seq(0, 1, 0.01)

for(is.reversible in c(0)){
  modelname <- "M9_IDLE"
  
  source(here("code", "initialize_model.R"))
  n_conditions <- 1
  rho_cond <- rho_cond[1]
  
  q0_wt <- q0
  last_feasible_q0 <- q0
  
  results_list <- list()
  idle <- which(reaction == "IDLE")
      
    for (fraction in phis_to_test){
      min_phi[idle] <- fraction
      max_phi[idle] <- fraction+1e-6
      
      source(here("code", "solver_loop.R"))
      
      fs <- q_opt[1, ]
      results_list[[length(results_list) + 1]] <- data.frame(
        phi = fraction,
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
  
  results <- do.call(rbind, results_list)
  write.csv(results, here("data", paste0(modelname, ".csv")))
}




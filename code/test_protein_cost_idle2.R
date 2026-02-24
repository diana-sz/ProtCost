# Test growth rate with two idle proteins expressed at once

rm(list = ls(all = TRUE))

library(here)

predict.parameters <- 0
phis_to_test <- seq(0.1, 0.4, 0.1)
idle1_fractions <- seq(0, 1, 0.2)

for (is.reversible in c(1, 0)) {
  modelname <- "M10_IDLE2"
  
  source(here("code", "initialize_model.R"))
  
  # Only one condition
  n_conditions <- 1
  rho_cond <- rho_cond[1]
  
  q0_wt <- q0
  last_feasible_q0 <- q0
  
  # Indices
  idle1_ind <- which(reaction == "IDLE")
  idle2_ind <- which(reaction == "IDLE2")
  
  results_list <- list()
  for(fraction in phis_to_test){
    for(f_idle1 in idle1_fractions){
      min_phi[idle1_ind] <- f_idle1 * fraction
      max_phi[idle1_ind] <- f_idle1 * fraction + 1e-6
      min_phi[idle2_ind] <- (1 - f_idle1) * fraction
      max_phi[idle2_ind] <- (1 - f_idle1) * fraction + 1e-6
      
      source(here("code", "solver_loop.R"))
      
      fs <- q_opt[1, ]
      results_list[[length(results_list) + 1]] <- data.frame(
        idle1_phi = f_idle1 * fraction,
        idle2_phi = (1 - f_idle1) * fraction,
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
  write.csv(results, here("data", paste0(modelname, ".csv")), row.names = FALSE)
}

# Test cost of an idle protein

rm(list=ls(all=TRUE))

library(here)

predict.parameters <- 0
#extreme_k <- 2000 # 2000  # FALSE
phis_to_test <- seq(0, 0.45, 0.01)

for(extreme_k in c(2000, FALSE)){
  for(is.reversible in c(0,1)){
    modelname <- "M10_Q_IDLE"
    
    source(here("code", "initialize_model.R"))
    n_conditions <- 1
    rho_cond <- rho_cond[1]
  
    q0_wt <- q0
    last_feasible_q0 <- q0
  
    if(extreme_k) KS[1,1] <- extreme_k
    
    results_list <- list()
    idle <- which(reaction == "IDLE")
    idle2 <- which(reaction == "Q")
    kcatf <- kcatf
    
      for (fraction in phis_to_test){
        print(paste("Testing fraction", fraction))
        min_phi[idle] <- fraction
        max_phi[idle] <- fraction+1e-6
        min_phi[idle2] <- 0.55
        
        source(here("code", "solver_loop.R"))
        
        fs <- q_opt[1, ]
        results_list[[length(results_list) + 1]] <- data.frame(
          phi = fraction,
          mu = mu_opt,
          convergence = res$convergence,
          t(c(
            setNames(fs, paste0("f.", reaction)),
            s          setNames(prot(fs), paste0("p.", reaction)),
            setNames(c(a, ci(fs)), paste0("c.", reactant))
          ))
        )
      }
    
    results <- do.call(rbind, results_list)
    
    suffix2 <- ifelse(extreme_k, "_extreme", "")
    write.csv(results, here("data", paste0(modelname, suffix2, ".csv")))
  }
}


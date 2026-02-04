# test the effect of toxic fuel production on growth rate and the expression 
# of an efflux pump

rm(list=ls(all=TRUE))

library(here)

directory <- paste0(here(), "/code")
setwd(directory) 
predict.parameters <- 0
phis_to_test <- seq(0, 0.5, 0.001)#, seq(0.5, 0, -0.005))
efflux_kcats <- c(0.1, 1, 10, 100)

for(is.reversible in c(1,0)){
  
  modelname <- "M10fuel_efflux"
  
  source("initialize_model.R")
  n_conditions <- 1
  
  q0_wt <- q0
  last_feasible_q0 <- q0
  
  results_list <- list()
  
  fuel <- which(reaction == "FUEL")
  efflux <- which(reaction == "EFFLUX")
  
  for(efflux_kcat in efflux_kcats){
    kcatf[efflux] <- efflux_kcat
    kcatb[efflux] <- is.reversible*efflux_kcat/5
    
    for (fuel_phi in phis_to_test){
      min_phi[fuel] <- fuel_phi
      max_phi[fuel] <- fuel_phi+1e-6
      
      source("solver_loop.R")
      
      fs <- q_opt[1, ]
      results_list[[length(results_list) + 1]] <- data.frame(
        kcat = efflux_kcat,
        fuel_phi = fuel_phi,
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
  
  results <- do.call(rbind, results_list)
  write.csv(results, paste0("../data/", modelname, "_tradeoff.csv"))
}

# Test the effect of secondary carbon source on the cost of the protein metabolizing it
# in the absence of x_L - LAC proteins have only cost, when x_L increases,
# the benefit of LAC also increases

rm(list=ls(all=TRUE))

library(here)

phis_to_test <- seq(0, 0.5, 0.005) # 0.005))
predict.parameters <- 0
primary_c_source <- 0.25
modelname_orig <- "M10_dekel_efflux" # "M10_dekel_efflux" M9_dekel

alt_concentrations <- c(0.001, 0.01, 0.05, 0.1, 0.2, 1, 10)

for(is.reversible in c(0)){
  results_list <- list()
  
  for (x_C2 in alt_concentrations){
    modelname <- modelname_orig
    source(here("code", "initialize_model.R"))
    q0_wt <- q0
    last_feasible_q0 <- q0
    
    alt_ind <- which(reaction == "LAC")
    a_cond[1,1] <- primary_c_source
    a_cond[3,1] <- x_C2 # alternative carbon source
    n_conditions <- 1

    source(here("code", "GBA_solver.R"))
    
    mu_orig <- mu_opt
    p_opt <- prot(q_opt[1,])
    
    for (fraction in phis_to_test){
      print(paste("LAC=", x_C2, "fraction=", fraction))
      
      min_phi[alt_ind] <- fraction
      max_phi[alt_ind] <- fraction+1e-6
      
      source(here("code", "solver_loop.R"))
      
      fs <- q_opt[1, ]
      results_list[[length(results_list) + 1]] <- data.frame(
        x_C = a_cond[1,1],
        x_C2 = x_C2,
        protein = "LAC",
        phi = fraction,
        mu = mu_opt,
        mu_norm = mu_opt / mu_orig,
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
  write.csv(results, here("data", paste0(modelname, ".csv")), row.names = FALSE)
  
}



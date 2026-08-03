# reproduce 'reverse growth laws' - the effect of ribosome inhibitor and 
# substrate concentration on growth rate

rm(list=ls(all=TRUE))

library(here)

predict.parameters <- 0
modelname_orig <- "M8_inh" 
inh_conc <- c(0.00001, 10^seq(-0.55, 1, by = 0.22))

for(is.reversible in c(0,1)){
  s_conc <- 10^seq(0.7, -1.1, by = -0.3)
  #if(is.reversible){
    s_conc <- 10^seq(0.7, -1, by = -0.25)
  #}
  
  
  results_list <- list()
  
  for (x_c in s_conc){
    
    for (x_inh in inh_conc){
        modelname <- modelname_orig
        source(here("code", "initialize_model.R"))
        q0_wt <- q0
        last_feasible_q0 <- q0
  
        a_cond[1,1] <- x_c
        a_cond[3,1] <- x_inh
        n_conditions <- 1
    
        source(here("code", "solver_loop.R"))
        
        fs <- q_opt[1, ]
        results_list[[length(results_list) + 1]] <- data.frame(
          x_C = x_c,
          x_INH = x_inh,
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
  write.csv(results, here("data", paste0(modelname, ".csv")))
}


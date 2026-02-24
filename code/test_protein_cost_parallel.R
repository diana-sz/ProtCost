# Test protein cost of each protein in the model
# run in parallel (one core - one reaction)

rm(list=ls(all=TRUE))

library(parallel)
library(here)


# --- read model and get initial solution --- ####
is.reversible <- 0
predict.parameters <- 0
modelname <- "M9_Q" # "M8" "M9_Q" "B"
max_cores <- 9

# Solver function ####
run_solver <- function(q0_value, msg = NULL) {
  if (!is.null(msg)) print(msg)
  assign("q0", q0_value, envir = .GlobalEnv)
  res_try <- try(source(here("code", "GBA_solver.R")), silent = TRUE)
  inherits(res_try, "try-error")
}


source(here("code", "initialize_model.R"))
n_conditions <- 1  # take only the first condition in the .ods model file
q0_wt <- q0

source(here("code", "GBA_solver.R"))

# --- prepare parameters for phi testing --- ####
results_list <- list()
p_opt <- prot(q_opt[1, ])
opt_phis <- (p_opt / rho) / sum(p_opt / rho)
mu_orig <- mu_opt

# --- Parallel setup --- ####
num_cores <- min(max_cores, detectCores() - 10)
batch_size <- ifelse(length(reaction) >= num_cores, num_cores, length(reaction))
reaction_batches <- split(reaction, ceiling(seq_along(reaction)/batch_size))

process_reaction <- function(reaction_name, opt_phi, q0_initial) {
  local_results <- list()
  last_feasible_q0 <- q0_initial
  last_q0 <- q0_initial
  
  
  # step size below optimum higher than above (the curves are steeper)
  stepsize <- 0.0001
  if(grepl("B", modelname)){
    stepsize <- 0.005
  }
  below <- seq(round(opt_phi, 4), 0, by = -stepsize)
  if(opt_phi/stepsize < 1000){
    below <- seq(round(opt_phi, 4), 0, length.out = 1000)

  }
  above <- seq(round(opt_phi, 4), 1, by = stepsize)
  
  phis_to_test <- c(below, above)

  for (fraction in phis_to_test) {
    print(paste(reaction_name, "phi =", fraction))
    
    # make local copies of bounds
    local_min_phi <- min_phi
    local_max_phi <- max_phi

    local_min_phi[reaction == reaction_name] <- fraction
    local_max_phi[reaction == reaction_name] <- fraction + 1e-6
    
    # assign to global env for solver
    assign("min_phi", local_min_phi, envir = .GlobalEnv)
    assign("max_phi", local_max_phi, envir = .GlobalEnv)
    
    # use last feasible q0 (for first simulation this is q0_wt)
    error <- FALSE
    error <- run_solver(
      last_feasible_q0,
      "Solving with last feasible solution"
    )
    
    # draw Gaussian noise with mean 0, sd = sd_vec
    # noise  <- rnorm(length(last_feasible_q0), mean = 0, sd = last_feasible_q0*0.02)
    # perturbed_q0 <- last_feasible_q0 + noise
    # perturbed_q0[perturbed_q0 < 0] <- 0

    candidates <- list(
      list(q = q0_wt,  msg = "Solver did not converge - trying with initial FBA solution"),
      list(q = q0_alt, msg = "Solver did not converge - trying with alternative FBA solution")
      #list(q = last_q0, msg = "Solver did not converge - trying with last solution")
    )
    
    for (cand in candidates) {
      if (res$convergence == -1){
        error <- run_solver(cand$q, cand$msg)
        if(error) print("Solver error")
      }
    }
    
    # if there is a converged solution, save the latest f0
    last_q0 <- q_initial
    if(res$convergence == 4){
      last_feasible_q0 <- q_initial
    }
    
    # do not save data if it does not converge to reduce final file size
    if(res$convergence != 4){
      next
    }
    
    qs <- q_opt[1, ]
    opt_phi <- ifelse(opt_phi < 1e-8, 1e-8, opt_phi)
    
    local_results[[length(local_results) + 1]] <- data.frame(
      reaction = reaction_name,
      phi = fraction,
      rel_phi = fraction/opt_phi,
      mu = mu_opt,
      mu_norm = mu_opt / mu_orig,
      convergence = res$convergence,
      t(c(
        setNames(qs, paste0("f.", reaction)),
        setNames(v(qs), paste0("v.", reaction)),
        setNames(prot(qs), paste0("p.", reaction)),
        setNames(c(a, ci(qs)), paste0("c.", reactant))
      ))
    )
  }
  
  do.call(rbind, local_results)
}

# --- Main loop for testing protein fractions --- ####
for (batch in reaction_batches) {
  cat("Running batch of size", length(batch), "...\n")
  
  results_batch <- mclapply(
    seq_along(batch),
    function(i) process_reaction(batch[i], opt_phis[reaction == batch[i]], q0),
    mc.cores = num_cores
  )
  
  results_list <- c(results_list, results_batch)
}

results <- do.call(rbind, results_list)
write.csv(results, here("data", paste0(modelname, "_protein_cost.csv")), row.names = FALSE)



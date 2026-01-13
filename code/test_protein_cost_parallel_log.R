rm(list=ls(all=TRUE))

library(parallel)
library(here)
library(rstudioapi)
library(readODS)
library(nloptr)
library(lpSolve)
library(Matrix)

set.seed(1)

directory <- file.path(here(), "code")
setwd(directory)

# Solver function ####
run_solver <- function(q0_value, msg = NULL) {
  if (!is.null(msg)) print(msg)
  assign("q0", q0_value, envir = .GlobalEnv)
  res_try <- try(source("GBA_solver_log.R"), silent = TRUE)
  inherits(res_try, "try-error")
}


# --- read model and get initial solution --- ####
is.reversible <- 1
predict.parameters <- 0
modelname <- "M8"
max_cores <- 8

source("initialize_model.R")
n_conditions <- 1
q0_wt <- q0

source("GBA_solver.R")

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
  
  #print(reaction_name, opt_phi)
  # opt_phi_nonzero <- TRUE
  # if (opt_phi < 1e-8) {
  #   opt_phi <- 0.1
  #   opt_phi_nonzero <- FALSE
  # }

  # step size below optimum higher than above (the curves are steeper)
  # go from both directions in case solver does not converge
  below <- seq(round(opt_phi, 4), 0, by = -0.0001) #ceiling(opt_phi/0.001))
  
  #problematic <- c("ATPS", "DNAP", "ADPS")
  #interval <- ifelse(reaction_name %in% problematic, 0.0005, 0.005)
  above <- seq(round(opt_phi, 4), 1, by = 0.0001)

  phis_to_test <- c(below, above)
  
  for (fraction in phis_to_test) {
    print(paste(reaction_name, "phi =", fraction))
    
    # make local copies of bounds
    local_min_phi <- min_phi
    local_max_phi <- max_phi
    
    if(fraction*opt_phi >= 1){next}
    local_min_phi[reaction == reaction_name] <- fraction
    local_max_phi[reaction == reaction_name] <- fraction + 1e-4
    
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
    noise  <- rnorm(length(last_feasible_q0), mean = 0, sd = last_feasible_q0*0.05)
    perturbed_q0 <- last_feasible_q0 + noise
    perturbed_q0[perturbed_q0 < 0] <- 0
      
    candidates <- list(
      list(q = perturbed_q0, msg = "Solver did not converge - trying with perturbed last solution"),
      list(q = last_q0, msg = "Solver did not converge - trying with last solution"),
      list(q = q0_wt,  msg = "Solver did not converge - trying with initial FBA solution"),
      list(q = q0_alt, msg = "Solver did not converge - trying with alternative FBA solution")
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
    
    qs <- q_opt[1, ]
    
    local_results[[length(local_results) + 1]] <- data.frame(
      x_Glc = a_cond[1, 1],
      reaction = reaction_name,
      phi = fraction,
      rel_phi = fraction/opt_phi,
      mu = mu_opt,
      mu_norm = mu_opt / mu_orig,
      #opt_phi_nonzero = opt_phi_nonzero,
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

write.csv(results, paste0("../data/", modelname, "_protein_cost_log.csv"))



rm(list=ls(all=TRUE))

library(parallel)
library(here)
library(rstudioapi)
library(readODS)
library(nloptr)
library(lpSolve)
library(Matrix)

directory <- file.path(here(), "code")
setwd(directory)

# --- read model and get initial solution --- ####
is.reversible <- 1
predict.parameters <- 0
modelname <- "A7simple"
max_cores <- 7

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
  #local_q0 <- q0_initial  # start with provided f0 for this reaction
  last_feasible_q0 <- q0_initial
  
  #print(reaction_name, opt_phi)
  opt_phi_nonzero <- TRUE
  if (opt_phi < 1e-8) {
    opt_phi <- 0.1
    opt_phi_nonzero <- FALSE
  }

  #rel_phis_to_test <- c(seq(1, 0.01, -0.01), seq(1, 1/opt_phi, 0.01))
  phis_to_test <- c(seq(0, 0.8, 0.002), seq(0.8, 0.0, -0.002))
  
  
  for (fraction in phis_to_test) {
    print(paste(reaction_name, "phi =", fraction))
    
    # make local copies of bounds
    local_min_phi <- min_phi
    local_max_phi <- max_phi
    
    if(fraction*opt_phi >= 1){next}
    local_min_phi[reaction == reaction_name] <- fraction
    local_max_phi[reaction == reaction_name] <- fraction + 1e-5
    
    # assign to global env for solver
    assign("min_phi", local_min_phi, envir = .GlobalEnv)
    assign("max_phi", local_max_phi, envir = .GlobalEnv)
    
    # if (fraction == 1){
    #   assign("q0", q0_wt, envir = .GlobalEnv)
    #   print("using wildtype q0")
    # }else{
    #   assign("q0", last_feasible_q0, envir = .GlobalEnv)
    #   print("using previous q0")
    # }
    
    # use last feasible q0 (for first simulation this is q0_wt)
    assign("q0", last_feasible_q0, envir = .GlobalEnv)
  
    error <- FALSE
    
    # try with the last solution
    error_check <- try({source("GBA_solver.R")}, silent = TRUE)
    if(class(error_check) == "try-error"){
      print("Solver error")
      error <- TRUE
    }
    
    # if it does not converge, try with the initial solution
    if(res$convergence == -1 | error){
      print("Solver did not converge - trying again with initial solution")
      # try again with the initial solution
      assign("q0", q0_wt, envir = .GlobalEnv)
      error_check <- try({source("GBA_solver.R")}, silent = TRUE)
      if(class(error_check) == "try-error"){
        print("Solver error")
        error <- TRUE
      }
    }
    
    # if it does not converge, try with alternative initial solution
    if(res$convergence == -1 | error){
      print("Solver did not converge - trying again with alternative initial solution")
      # try again with the initial solution
      assign("q0", q0_alt, envir = .GlobalEnv)
      error_check <- try({source("GBA_solver.R")}, silent = TRUE)
      if(class(error_check) == "try-error"){
        print("Solver error")
        error <- TRUE
      }
    }
    
    # if there is a converged solution, save the latest f0
    if(res$convergence == 4){
      last_feasible_q0 <- q_initial
    }
    
# 
#     # Always update local_f0 regardless of convergence
#     local_q0 <- q0  
    
    qs <- q_opt[1, ]
    
    local_results[[length(local_results) + 1]] <- data.frame(
      x_Glc = a_cond[1, 1],
      reaction = reaction_name,
      phi = fraction,
      rel_phi = fraction/opt_phi,
      mu = mu_opt,
      mu_norm = mu_opt / mu_orig,
      opt_phi_nonzero = opt_phi_nonzero,
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

write.csv(results, paste0("../data/", modelname, is.reversible, "_protein_cost.csv"))


# --- Plot results and save as pdf --- ####
# 
# MAX_Y <- 1.1
# 
# pdf(
#   file.path("../figures", paste0(modelname, "_conditions_cost_test_", is.reversible, ".pdf")),
#   title = paste("Protein cost testing", modelname),
#   width = 10, height = 5
# )
# 
# par(mar = c(5, 5, 5, 1), mfrow = c(1, 2))
# 
# for (rxn in unique(results$reaction)) {
#   one_rxn <- subset(results, reaction == rxn)
#   colors <- as.factor(one_rxn$x_Glc)
#   
#   # Plot normalized growth
#   with(subset(one_rxn, opt_phi_nonzero == TRUE),
#        plot(mu_norm ~ rel_phi,
#             xlab = "Titrated level / opt. level",
#             ylab = "Normalized growth (mu)",
#             ylim = c(0, MAX_Y), xlim = c(0, 10),
#             pch = shape, cex = 1.3, cex.lab = 1.3,
#             main = rxn, col = colors))
#   
#   # Plot absolute growth
#   with(one_rxn,
#        plot(mu ~ phi,
#             xlab = "Proteome fraction",
#             ylab = expression("Growth rate" ~ "[" * h^-1 * "]"),
#             ylim = c(0, max(one_rxn$mu)*1.1), xlim = c(0, max(one_rxn$phi)*1.1),
#             pch = shape, cex = 1.3, cex.lab = 1.3,
#             main = rxn, col = colors))
# }
# 
# dev.off()

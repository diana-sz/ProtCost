
error <- FALSE

# try with the last solution
q0 <- last_feasible_q0
error_check <- try({source("GBA_solver.R")}, silent = TRUE)
if(class(error_check) == "try-error"){
  print("Solver error")
  error <- TRUE
}

# if it does not converge, try with the initial solution
if(res$convergence == -1 | error){
  print("Solver did not converge - trying again with initial solution")
  # try again with the initial solution
  q0 <- q0_wt
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
  q0 <- q0_alt
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
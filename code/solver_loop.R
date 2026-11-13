# nloptr codes: 1 to 4 = success, 5 and 6 = limit reached, negative = failure

# try with the last solution
q0 <- last_feasible_q0
error_check <- try({source(here("code", "GBA_solver.R"))}, silent = TRUE)
error <- inherits(error_check, "try-error")
if (error) print("Solver error")

# if it does not converge, try with the initial solution
if (error || res$convergence < 0) {
  print("Solver did not converge - trying again with initial solution")
  q0 <- q0_wt
  error_check <- try({source(here("code", "GBA_solver.R"))}, silent = TRUE)
  error <- inherits(error_check, "try-error")
  if (error) print("Solver error")
}

# if it does not converge, try with alternative initial solution
if (error || res$convergence < 0) {
  print("Solver did not converge - trying again with alternative initial solution")
  q0 <- q0_alt
  error_check <- try({source(here("code", "GBA_solver.R"))}, silent = TRUE)
  error <- inherits(error_check, "try-error")
  if (error) print("Solver error")
}

# if there is a converged solution, save the latest q0
if (!error && res$convergence == 4) {
  last_feasible_q0 <- q_initial
}
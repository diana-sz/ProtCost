# q0 ########################################################
# estimates q0 by minimizing total protein at given mu_data 
a <- a_cond[,1]
rho <- rho_cond[1]

initial_q <- function(lowerp) {

  lower_met <- biomass_lb[-p]

  lowerb <- c(lower_met, lowerp)
  
  q0 <- rep(0,r)
  
  # defining parameters
  objective.fn <- rep(1,r)   # to minimize the sum q/kcat
  
  const.mat <- rbind(M,s,diag(r)) # flux balance, density, non-negative q
  
  const.dir <- c(rep(">=",p),"=",rep(">=",r) )
  
  const.rhs <- c(lowerb, 1, min_q)
  
  # solving model
  
  lp.solution <- lp("min", objective.fn, const.mat, const.dir, const.rhs, timeout = 10)
  
  # solution to the linear optimization
  
  q0 <- lp.solution$solution
  
  n_tests <- 0
  while (sum(q0) == 0 & n_tests < 20) {
    n_tests <- n_tests + 1
    
    lowerb[-length(lowerb)] <- lowerb[-length(lowerb)]/1.1
    
    const.rhs <- c(lowerb, 1, min_q)
    
    # solving model
    lp.solution <- lp("min", objective.fn, const.mat, const.dir, const.rhs, timeout = 10)
    
    # solution to the linear optimization
    q0 <- lp.solution$solution
    
  }
  
  return(q0)
}
  

q0 <- initial_q(biomass_lb[p])

# finds best q0 by testing different bp_min
if (best_q0 == T) {
  
  # mu ############## (growth rate)
  mu <- function(q) as.numeric( M[p,r]*q[r]/(tau(a,rho*b(q))%*%q)  )
  
  # biomass fractions
  b <- function(q) M%*%q
  
  n_tests <- biomass_lb[p]/0.005
  
  # find best q0 by varying cp_min
  q0_test <- matrix(rep(0,r*n_tests),nrow = n_tests)
  for (test in 1:n_tests) q0_test[test,] <- initial_q(biomass_lb[p]*0.5 + 0.005*test)
  
  mu_cp <- rep(0,n_tests)
  for (test in 1:n_tests) mu_cp[test] <- mu(q0_test[test,])
  
  q0_test <- q0_test[!is.nan(mu_cp),]
  
  mu_cp   <- mu_cp[!is.nan(mu_cp)]

  q0_best <- q0_test[mu_cp==max(mu_cp),]
  
  if (sum(q0_best) != 0) q0 <- q0_best

}




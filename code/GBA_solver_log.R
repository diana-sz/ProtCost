library(nloptr)

# Functions ####################################################################

# mu ############## (growth rate)
mu <- function(q) as.numeric( (M[p,r]*q[r]/(tau(a,rho*b(q))%*%q) ) )

# negative log mu (for minimization)
negative_log_mu <- function(q) {
  m <- mu(q)
  if (m <= 0 || !is.finite(m)) return(1e12)
  -log(m)
}

# fluxes
v <- function(q) as.vector(mu(q))*rho*q

# protein concentrations
prot <- function(q) tau(a,rho*b(q))*v(q)

# internal concentrations "i" of metabolites "m" and proteome "p"
ci <- function(q) rho*M%*%q

# biomass fractions
b <- function(q) M%*%q

# proteome fractions
phi <- function(q) (tau(a,rho*b(q))*v(q))/(ci(q)[p])

# Optimization #################################################################

# equality constraints (density constraint)
g <- function(q) s%*%q - 1

# equality constraints if also constraint on ribosome composition
if (ribcomp > 0 & "rRNA" %in% i_reactant)
  g <- function(q) c(s%*%q - 1, rprna(q) - ribcomp )

# inequality constraints (min c, max c, min phi, max phi)
h <- function(q) c( ci(q) - min_c, max_c - ci(q), phi(q) - min_phi, max_phi - phi(q))

# Indirect elasticities ########################################################

E <- function(q) rho*dtau(a,rho*b(q))%*%M

# Gradient of negative mu ######################################################

negative_dmu <- function(q) {
  m <- mu(q)
  if (m <= 0 || !is.finite(m)) return(rep(0, length(q)))
  -as.numeric( ((m^2)/b(q)[p]) *
                 (M[p,]/m -
                    t(q)%*%(rho*dtau(a,rho*b(q))%*%M) -
                    tau(a,rho*b(q))) )
}

# Gradient of negative log mu ##################################################

negative_dlogmu <- function(q) {
  m <- mu(q)
  if (m <= 0 || !is.finite(m)) return(rep(0, length(q)))
  negative_dmu(q) / m
}

# starts loop for the optimization at each environmental condition #############

q_opt  <- matrix(rep(0,r*n_conditions),ncol=r)
mu_opt <- rep(0,n_conditions)

otime  <- rep(0,n_conditions)
conv   <- rep(0,n_conditions)
iter   <- rep(0,n_conditions)

q_initial <- q0

# Upper and lower bounds
upper_q <- max_q
lower_q <- min_q

suppressMessages(
  for (cond in 1:n_conditions) {
    
    a <- a_cond[,cond]
    rho <- rho_cond[cond]
    
    st <- system.time({
      res <- auglag(
        q_initial,
        fn = negative_log_mu,
        gr = negative_dlogmu,
        lower = lower_q,
        upper = upper_q,
        heq = g,
        hin = h,
        localsolver = c(solver),
        localtol = mu_tol,
        control = list(
          maxeval = 100000,
          xtol_rel = xtol_rel
          #xtol_abs = xtol_abs
        )
      )
    })
    
    opt_time <- signif(st[[1]], digits = 4)
    
    # update initial condition
    q_initial <- res$par
    
    # store results
    q_opt[cond,] <- q_initial
    mu_opt[cond] <- mu(q_initial)
    conv[cond]   <- res$convergence
    otime[cond]  <- opt_time
    iter[cond]   <- res$iter
    
    print(paste(
      "optimization: ", cond, "/", n_conditions,
      ", time: ", otime[cond], "s",
      ", iterations: ", iter[cond],
      ", convergence: ", conv[cond],
      ", growth rate: ", as.numeric(signif(mu_opt[cond], digits = 5)),
      sep = ""
    ))
  }
)

# produces c_opt ###############################################################

c_opt <- matrix(rep(0,p*n_conditions),ncol=p)
for (cond in 1:n_conditions){
  rho <- rho_cond[cond]
  c_opt[cond,] <- as.numeric(t(ci(q_opt[cond,])))
}

# mean optimization time #######################################################

mean_time <- signif(mean(otime), digits = 3)

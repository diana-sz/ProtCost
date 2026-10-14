library(nloptr)

# Functions ####################################################################

# mu ############## (growth rate)
mu <- function(q) as.numeric( (M[p,r]*q[r]/(tau(a,rho*b(q))%*%q) ) )

# negative mu (for minimization)
negative_mu <- function(q) -as.numeric(M[p,r]*q[r]/(tau(a,rho*b(q))%*%q)) 

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


# Density as a function of growth rate #########################################

# rho = rho_0 - rho_slope * mu
# rho_0 is the density at mu = 0, not the density of condition 1
rho_0     <- rho_cond[1]
rho_slope <- 100
rho_fun   <- function(mu) rho_0 - rho_slope * mu

# fixed-point loop for rho
omega        <- 0.5    # damping factor, lower it if rho and mu oscillate
rho_rel_tol  <- 1e-6   # relative tolerance for rho
rho_max_iter <- 50     # maximum number of rho updates per condition


# Optimization #################################################################

# equality constraints (density constraint)
g <- function(q) s%*%q - 1

# equality constraints if also constraint on ribosome composition
if (ribcomp > 0 & "rRNA" %in% i_reactant) g <- function(q) c(s%*%q - 1, rprna(q) - ribcomp )

# inequality constraints (min c, max c, min phi, max phi)
h <- function(q) c( ci(q) - min_c, max_c - ci(q), phi(q) - min_phi, max_phi - phi(q))

# Indirect elasticities ########################################################

E <- function(q) rho*dtau(a,rho*b(q))%*%M

# Gradient of negative_mu ######################################################

negative_dmu <- function(q) -as.numeric( ((mu(q)^2)/b(q)[p])*(M[p,]/mu(q) - 
                                                                t(q)%*%(rho*dtau(a,rho*b(q))%*%M) - tau(a,rho*b(q))) )

# starts loop for the optimization at each environmental condition #############
q_opt  <- matrix(rep(0,r*n_conditions),ncol=r) 
mu_opt <- rep(0,n_conditions) 

otime    <- rep(0,n_conditions)
conv     <- rep(0,n_conditions)
iter     <- rep(0,n_conditions)
rho_opt  <- rep(0,n_conditions)
rho_iter <- rep(0,n_conditions)
dmu_opt  <- matrix(rep(0,r*n_conditions),ncol=r)

q_initial <- q0

# Upper bounds 
upper_q <- max_q

# lower bounds 
lower_q <- min_q


suppressMessages(
  for (cond in 1:n_conditions) {
    
    a      <- a_cond[,cond]
    rho    <- rho_cond[cond]   # starting value for rho
    rho_ok <- FALSE
    
    st <- system.time({
      
      for (k in 1:rho_max_iter) {
        
        # optimization with fixed rho
        res <- auglag(q_initial, fn = negative_mu, lower = lower_q, gr = negative_dmu,
                      upper = upper_q, heq = g, hin = h,
                      localsolver = c(solver), localtol = mu_tol,
                      control = list(maxeval = 100000, xtol_rel = xtol_rel))
        
        # warm start for the next solve
        q_initial <- res$par
        
        # density implied by the current growth rate
        rho_new <- rho_fun(mu(q_initial))
        if (rho_new <= 0) stop(paste("rho <= 0 in condition", cond,
                                     ": mu is too high for this rho_fun"))
        
        # stop if rho is consistent with mu
        if (abs(rho_new - rho) <= rho_rel_tol * rho) { rho_ok <- TRUE; break }
        
        # damped update
        rho <- rho + omega * (rho_new - rho)
      }
      
    })
    
    if (!rho_ok) warning(paste("rho did not converge in condition", cond))
    
    # convergence (codes: "-1" means optimization problem, "5" means optimization 
    # stop because maxeval, "4" means optimization stopped because xtol_rel or 
    # xtol_abs (above) was reached)
    
    opt_time  <- signif(st[[1]], digits= 4)
    
    # solution
    q_opt[cond,] <- q_initial
    
    # converged density
    rho_opt[cond]  <- rho
    rho_iter[cond] <- k
    
    # optimal mu
    mu_opt[cond] <- mu(q_initial)
    
    conv[cond] <- res$convergence
    
    # optimization time
    otime[cond] <- opt_time
    
    # number of iterations (last solve)
    iter[cond] <- res$iter  
    
    print(paste("optimization: ",cond, "/", n_conditions,", time: ",
                otime[cond],"s, iterations: ",iter[cond],", rho updates: ",rho_iter[cond],
                ", convergence: ", conv[cond],", rho: ", signif(rho_opt[cond], digits= 5),
                ", growth rate: ", as.numeric(signif(mu_opt[cond], digits= 5)), sep=""))
  }
)


# produces c_opt (uses the converged rho of each condition)
c_opt <- matrix(rep(0,p*n_conditions),ncol=p)
for (cond in 1:n_conditions){
  a   <- a_cond[,cond]
  rho <- rho_opt[cond]
  c_opt[cond,] <- as.numeric(t(ci(q_opt[cond,])))
}  

# check for active bounds in h (values close to 0 mean an active bound)
active <- NULL

for (cond in 1:n_conditions){
  a   <- a_cond[,cond]
  rho <- rho_opt[cond]
  q   <- q_opt[cond,]
  
  hv    <- h(q)
  n_c   <- length(ci(q))
  n_phi <- length(phi(q))
  
  block <- rep(c("c - min_c", "max_c - c", "phi - min_phi", "max_phi - phi"),
               times = c(n_c, n_c, n_phi, n_phi))
  
  nm_c   <- if (is.null(rownames(M))) seq_len(n_c)   else rownames(M)
  nm_phi <- if (is.null(colnames(M))) seq_len(n_phi) else colnames(M)
  name   <- c(nm_c, nm_c, nm_phi, nm_phi)
  
  idx <- which(hv < 1e-6)
  if (length(idx) > 0)
    active <- rbind(active, data.frame(cond = cond, index = idx, block = block[idx],
                                       name = name[idx], value = hv[idx]))
}

print(active)

# mean optimization time
mean_time <- signif(mean(otime), digits= 3)
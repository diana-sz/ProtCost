library(Matrix)
library(lpSolve)

modelname <<- modelname

# Predicts kinetic parameters
predict.parameters <<- predict.parameters

# Forces all reactions to be irreversible if desired
is.reversible <<- is.reversible

# verbose: prints optimization diagnostics
verbose <<- F

# tolerances for the solver
xtol_rel <<- 1e-10
mu_tol <<- 1e-13

# solver 
solver <<- "SLSQP" # "SLSQP" # "gradi" LBFGS


suppressMessages(source(here("code", "Readmodelods.R")))

if (is.reversible == 1) modelname <- paste(modelname,"_rev",sep="")


# kinetics #####################################################################
source(here("code", "Kinetics.R"))


# Finds initial condition ##########################################
# finds best q0 first by testing different b_p
best_q0 <<- F
source(here("code", "q0_biomass.R"))
q0_alt <- q0

if (predict.parameters == 0){
  best_q0 <<- T
  source(here("code", "q0_biomass.R"))
}


# Predicts kinetic parameters based on mu and phi data #########################

if (predict.parameters > 0 & sum(q0) != 0) source(here("code", "Parameter_prediction.R"))


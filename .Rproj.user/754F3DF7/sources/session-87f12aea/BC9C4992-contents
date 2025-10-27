# Interface

# Clear variables
rm(list=ls(all=TRUE))

require('rstudioapi') 
require('readODS')
require('nloptr')
require('Matrix')
require('MASS')
require('lpSolve')
require("latex2exp")
require("lattice")
require("deSolve")
require("rgl")
require("plot3D")
#require("car")
require('gplots')
require('Rlibeemd')
#require('sfsmisc')
#require('plotly')

# Sets working directory to source file location ###############################

directory <- dirname(getActiveDocumentContext()$path)

setwd(directory) 

# Growth Balance Analysis (GBA) function for balanced growth simulations #######

GBA <- function(modelname,predict.parameters,is.reversible) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible")))
  
  modelname <<- modelname
  
  # Predicts kinetic parameters
  predict.parameters <<- predict.parameters
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  source("GBA.R")
  
}

# Example:
# GBA("B13",2,1)

# Growth Mechanics (GM) function for dynamic simulations #######################

GM <- function(modelname,predict.parameters,is.reversible,delta_t,totalT) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible",
                  "delta_t","totalT")))
  
  modelname <<- modelname
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  # Predicts kinetic parameters if > 0, using cm/Km = predict.parameters
  predict.parameters <<- predict.parameters
  
  # delta t
  delta_t <<- delta_t
  
  # Total simulation time
  totalT <<- totalT
  
  # reads model
  suppressMessages(source("Readmodelods.R"))
  
  if (is.reversible == 0) modelname <<- paste(modelname,"i",sep="")
  
  rtol <<- 1e-8 
  
  # defines dummy jo = 0, this variable is only accounted for in the other function
  # GMio
  
  jo <<- 0
  
  suppressMessages(source("GM.R"))
  
}

# Example:
#GM("L3per",3,0,0.001,3)
  
# Function to do many GM simulations with different x, given a model with fixed x
GMmeta <- function(modelname,predict.parameters,is.reversible,delta_t,totalT) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible",
                          "delta_t","totalT")))
  
  modelname <<- modelname
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  # Predicts kinetic parameters
  predict.parameters <<- predict.parameters
  
  # delta t
  delta_t <<- delta_t
  
  # Total simulation time
  totalT <<- totalT
  
  jo <<- 0
  
  # reads model
  suppressMessages(source("Readmodelods.R"))
  
  if (is.reversible == 0) modelname <<- paste(modelname,"i",sep="")
  
  # get x such that tau are equally spaced, better plots
  taumax <- (1/6)*(1 + 1/0.7)
    
  taumin <- 1/6 + 0.001
  
  tauseq <- seq(taumax,taumin,-(taumax-taumin)/5)
  
  xmeta <- 1/(6*tauseq -1)
    
  results_freq <- matrix(rep(0,length(xmeta)*5),ncol=5)
  
  for (itest in 1:length(xmeta)) {
    
    #x <<- 100/itest
    
    xt <<- function(t) xmeta[itest]
    
    rtol <<- 1e-8
    
    source("GM.R")
    
    # results to save (x,freq, mu, amplitude)
    results_freq[itest,] <- c(x, mean(mu_opt),mu_f,  max(mu_opt) - min(mu_opt),mu_E0)
    
    if(itest == 1) results_phir <- cbind(itest,mu_opt,phit[,r])
    
    if(itest > 1) results_phir  <- rbind(results_phir, cbind(itest,mu_opt,phit[,r]))
    
  }
    
    setwd(paste(directory,"/Meta",sep=""))
    
    # export results
    write.csv(results_freq, file = paste("frequencies ",modelname,".csv",sep=""))
    
    write.csv(results_phir, file = paste("phir ",modelname,".csv",sep=""))
    
    setwd(directory)
  
}

# Example:
# GMmeta("Aper",1,0,0.001,10)

results_freq <- list(0)
results_mag  <- list(0) 

# Function to do many GM simulations with different x, given a model with fixed x
GMresonance <- function(modelname,predict.parameters,is.reversible,delta_t,totalT) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible",
                          "delta_t","totalT")))
  
  modelname <<- modelname
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  # Predicts kinetic parameters
  predict.parameters <<- predict.parameters
  
  # delta t
  delta_t <<- delta_t
  
  # Total simulation time
  totalT <<- totalT
  
  # reads model
  suppressMessages(source("Readmodelods.R"))
  
  if (is.reversible == 0) modelname <<- paste(modelname,"i",sep="")
  
  nenv <- 5
  
  # xi_E environment frequency
  xi_E <<- 3*exp(-2:2)
 
  for (itest in 1:nenv) {
    
    xt   <<- function(t) 10 + 5*sin(xi_E[itest]*t) 
    
    dxdt <<- function(t) 5*xi_E[itest]*cos(xi_E[itest]*t) 
    
    rtol <<- 1e-8
    
    source("GM.R")
    
    # results to save (x,freq, mu, amplitude)
    results_freq[[itest]] <<- frequencies[-1]
    
    results_mag[[itest]]  <<- magnitudes[-1]
    
  }
  
}

# Example:
# GMresonance("Asin",3,0,0.001,100)

GMio <- function(modelname,predict.parameters,is.reversible,delta_t,totalT,jo,omega) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible",
                          "delta_t","totalT","jo","omega")))
  
  modelname <<- modelname
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  # Predicts kinetic parameters if > 0, using cm/Km = predict.parameters
  predict.parameters <<- predict.parameters
  
  # delta t
  delta_t <<- delta_t
  
  # Total simulation time
  totalT <<- totalT
  
  # jo = reaction j to have kcatj oscillated
  jo <<- jo
  
  # angular frequency of the imposed oscillation
  omega <<- omega
  
  # reads model
  suppressMessages(source("Readmodelods.R"))
  
  if (is.reversible == 0) modelname <<- paste(modelname,"i",sep="")
  
  rtol <<- 1e-8
  
  suppressMessages(source("GM.R"))
  
}

# Example:
# GMio("A",3,0,0.001,10,2,1/pi)



#dev.off()


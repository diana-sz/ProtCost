rm(list=ls(all=TRUE))

library(here)
library(xtable)

directory <- paste0(here(), "/code")
setwd(directory)

predict.parameters <- 0

for(is.reversible in c(1,0)){
  modelname <- "A7simple"

  source("initialize_model.R")
  
  source("GBA_solver.R") 
  
  source("GBA_Exportcsv.R")
  
  source("GBA_Plots_v2.R")
}





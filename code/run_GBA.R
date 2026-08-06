rm(list=ls(all=TRUE))

library(here)

predict.parameters <- 0

for(is.reversible in c(1,0)){
  modelname <- "M10_Q_IDLE"

  source(here("code", "initialize_model.R"))
  
  source(here("code", "GBA_solver.R"))
  
  source(here("code", "GBA_Exportcsv.R"))
}





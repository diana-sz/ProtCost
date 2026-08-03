# Script to predict kinetic parameters

rm(list=ls(all=TRUE))

library(here)
library(xtable)

modelname <- "M9_Q"
rescale_kcats <- TRUE
predict.parameters <- 3  # saturation parameter used in 'initialize_model.R'
is.reversible <- 1

source(here("code", "initialize_model.R"))

print("Predicted parameters:")
print(kcatf)
print(kcatb)
print(K)


latex_table <- xtable(rbind(kcatf, kcatb), 
                      caption = "Predicted turnover numbers", 
                      label = "tab:kcats",
                      digits = 1)
print(latex_table)


latex_table <- xtable(K, 
                      caption = "Predicted Km", 
                      label = "tab:kms",
                      digits = 1)
print(latex_table)


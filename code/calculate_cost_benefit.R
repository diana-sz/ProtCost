library(here)
library(readODS)
library(rstudioapi)
library(xtable)

directory <- paste0(here(), "/code")
setwd(directory)


modelname="A7simple"
suppressMessages(source("Readmodelods_v2.R"))

source("GBA_Kinetics.R")


opt_data <- read.csv("Results GBA/GBA Model A7simple mean time (5.11s) results.csv", row.names = 1)

row <- 1
rho <- rho_cond[1]
x  <- x_cond[,1]

vs <- opt_data[row, grep("v\\.", colnames(opt_data))]
fs <- opt_data[row, grep("f\\.", colnames(opt_data))]

cint <- opt_data[row, which(colnames(opt_data)=="C"):which(colnames(opt_data)=="p")]
fint <- cint/rho
taus <- tau(t(cint))
dtaus <- dtau(t(cint))
growth_rate <- opt_data[row, "mu"]

results <- data.frame(
  variable = character(),
  prot = numeric(),
  local_cost = numeric(),
  local_benefit = numeric(),
  transport_benefit = numeric(),
  sum = numeric(),
  stringsAsFactors = FALSE
)


for(j in 1:length(vs)){
  Mjp <- M["p",j]
  local_cost <- growth_rate*taus[j]
  local_benefit <- unlist(vs) %*% (dtaus %*% M[, j])
  transport_benefit <- colSums(M)[j] * unlist(vs) %*% dtaus %*% unlist(fint)
  

  print(paste(colnames(vs)[j],
              "prot", Mjp,
              "local cost", round(local_cost, 4),
              "local benefit", round(local_benefit, 4),
              "transport benefit", round(transport_benefit, 4),
              "sum", round(Mjp  - local_cost - local_benefit + transport_benefit, 4)))  
  
  results <- rbind(
    results,
    data.frame(
      variable = colnames(vs)[j],
      prot = as.numeric(Mjp),
      local_cost = -as.numeric(local_cost),
      local_benefit = -as.numeric(local_benefit),
      transport_benefit = as.numeric(transport_benefit),
      sum = as.numeric(Mjp - local_cost - local_benefit + transport_benefit)
    )
  )
}

latex_table <- xtable(results, 
                      caption = "GBA Model Results", 
                      label = "tab:GBA_results",
                      digits = 2)
print(latex_table)

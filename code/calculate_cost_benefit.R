# Calculate cost and benefit contributions of each reaction to the total growth
# cost at optimal growth rate
rm(list=ls(all=TRUE))

library(here)
library(xtable)

modelname <- "M8"
is.reversible <- 1

suppressMessages(source(here("code", "Readmodelods.R")))
source(here("code", "Kinetics.R"))

opt_data <- read.csv(here("data", "M8_rev_GBA.csv"), row.names = 1)

row <- 1
rho <- rho_cond[1]
a  <- a_cond[,1]

vs <- opt_data[row, grep("v\\.", colnames(opt_data))]
fs <- opt_data[row, grep("f\\.", colnames(opt_data))]

cint <- opt_data[row, which(colnames(opt_data)=="S"):which(colnames(opt_data)=="P")]
fint <- cint/rho
taus <- tau(a, t(cint))
dtaus <- dtau(a, t(cint))
growth_rate <- opt_data[row, "mu"]

results <- data.frame(
  variable = character(),
  prot = numeric(),
  catalytic_cost = numeric(),
  saturation_value = numeric(),
  crowding_value = numeric(),
  sum = numeric(),
  stringsAsFactors = FALSE
)


for(j in 1:length(vs)){
  Mjp <- M["P",j]
  catalytic_cost <- growth_rate*taus[j]
  saturation_value <- unlist(vs) %*% (dtaus %*% M[, j])
  crowding_value <- colSums(M)[j] * unlist(vs) %*% dtaus %*% unlist(fint)
  
  results <- rbind(
    results,
    data.frame(
      variable = colnames(vs)[j],
      prot = as.numeric(Mjp),
      catalytic_cost = -as.numeric(catalytic_cost),
      saturation_value = -as.numeric(saturation_value),
      crowding_value = as.numeric(crowding_value),
      sum = as.numeric(Mjp - catalytic_cost - saturation_value + crowding_value)
    )
  )
}

rownames(results) <- gsub("v.", "", results$variable,)

latex_table <- xtable(results[,2:ncol(results)], 
                      caption = "GBA Model Results", 
                      label = "tab:GBA_results",
                      digits = 3)
print(latex_table)

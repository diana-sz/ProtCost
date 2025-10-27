library(here)
library(readODS)
library(rstudioapi)
library(xtable)
#library(RColorBrewer)

directory <- paste0(here(), "/code")
setwd(directory)


modelname="A7simple"
is.reversible <- 1
predict.parameters <- 0

source("uni_colors.R")
source("initialize_model.R")

opt_data <- read.csv(paste0("../data/", modelname, is.reversible, "_protein_cost.csv"), row.names = 1)
opt_data <- opt_data[opt_data$convergence == 4, ]

row <- 1
rho_cond <- rho_cond[1]
a_cond <- a_cond[,1, drop=FALSE]


results <- data.frame(
  growth_rate = numeric(),
  tested_protein = character(),
  phi = numeric(),
  variable = character(),
  protein_benefit = numeric(),
  local_cost = numeric(),
  local_benefit = numeric(),
  transport_benefit = numeric(),
  sum = numeric(),
  stringsAsFactors = FALSE
)

for (tested_protein in unique(opt_data$reaction)) {
  one_prot <- opt_data[opt_data$reaction == tested_protein,]
  for (phi in one_prot$phi){
    one_phi <- one_prot[one_prot$phi == phi,]
    
    vs <- one_phi[row, grep("v\\.", colnames(one_phi))]
    fs <- one_phi[row, grep("f\\.", colnames(one_phi))]
    
    cint <- one_phi[row, grep("c\\.", colnames(one_phi))]
    cint <- cint[, (1+length(a_cond)):length(cint)]
    fint <- cint/rho
    taus <- tau(0, t(cint))
    dtaus <- dtau(0, t(cint))
    growth_rate <- one_phi[row, "mu"]
    
    for(j in 1:length(vs)){
      Mjp <- M["p",j]
      local_cost <- growth_rate*taus[j]
      local_benefit <- unlist(vs) %*% (dtaus %*% M[, j])
      transport_benefit <- colSums(M)[j] * unlist(vs) %*% dtaus %*% unlist(fint)
      
      results <- rbind(
        results,
        data.frame(
          growth_rate = growth_rate,
          tested_protein = tested_protein,
          phi = phi,
          variable = colnames(vs)[j],
          protein_benefit = as.numeric(Mjp),
          local_cost = -as.numeric(local_cost),
          local_benefit = -as.numeric(local_benefit),
          transport_benefit = as.numeric(transport_benefit),
          sum = as.numeric(Mjp - local_cost - local_benefit + transport_benefit)
        )
      )
    }
  }
}
  

xlim <- c(0,0.8)
ylim <- c(-2.1,2.1)
plots <- c("local_cost", "local_benefit", "transport_benefit", "sum")
#colors <- brewer.pal(length(plots), "Set1")
colors <- c(uni_red, uni_green, uni_lila, uni_blue)

n_prot <- length(unique(results$tested_protein))
n_rxn  <- length(unique(results$variable))

png(paste0("../figures/", modelname, is.reversible, "_cost_benefit.png"), 
    type="cairo", units="cm",
    width=22, height=22, res=300)

par(mfrow = c(n_rxn + 1, n_prot), mar = c(1.7,0,0,0.5), oma = c(5,5.5,0.5,0.5), xpd = FALSE)

# plot growth rates on top
prot_idx <- 0
for(protein in unique(results$tested_protein)){
  prot_idx <- prot_idx + 1
  one_protein <- results[results$tested_protein == protein, ]
  plot(one_protein$phi, one_protein$growth_rate, axes = FALSE, 
       pch = 19, cex = 0.9,
       ylim = c(0, max(one_protein$growth_rate)*1.05), xlim = xlim)
  abline(v=one_protein$phi[which.max(one_protein$growth_rate)], col="grey70", lty=2)
  if(prot_idx == 1){
    axis(2, at = c(0, 0.5, 1))
    mtext("Growth rate", side = 2, line = 2.5, cex = 1)
  }
  box()
}

# plot cost-benefit curves
par(mar = c(0,0,0.5,0.5))
rxn_idx <- 0
for(reaction in rev(unique(results$variable))){
  rxn_idx <- rxn_idx + 1
  one_rxn <- results[results$variable == reaction, ]
  
  prot_idx <- 0
  for(prot in unique(one_rxn$tested_protein)){
    prot_idx <- prot_idx + 1
    one_prot <- one_rxn[one_rxn$tested_protein == prot, ]
    
    # plot background
    plot(NA, xlim = xlim, ylim = ylim,
         axes = FALSE, xlab = "", ylab = "")
    abline(h=0)
    abline(v=one_prot$phi[which.max(one_prot$growth_rate)], col="grey70", lty=2)
    box()
    
    # add all curves
    for(plot_idx in seq_along(plots)){
      # points(one_prot$phi, one_prot[[plots[plot_idx]]],
      #        pch = 19, col = colors[plot_idx], cex=0.6)
      lines(one_prot$phi, one_prot[[plots[plot_idx]]],
             lwd = 2, col = colors[plot_idx])
    }
    
    # only y-axis on first column
    if(prot_idx == 1){
      axis(2, at = c(-2,-1,0,1,2))
      mtext(gsub("v\\.", "", reaction), side = 2, line = 2.5, cex = 0.85)
    }
    
    # only x-axis on last row
    if(rxn_idx == n_prot){
      axis(1)
      mtext(prot, side = 1, line = 2.5, cex = 0.85)
    }
    
    # add legend once in outer margin
    if(rxn_idx == 1 & prot_idx == 1){
      par(xpd = NA)  # allow legend outside plot region
      legend(0, 3.6, horiz = TRUE,
             legend = gsub("_", " ", plots), 
             col = colors, 
             pch = 19, bty = "n", cex = 1.2)
      par(xpd = FALSE)
      
    }

  }
}

# Add one common x and y label
mtext("Proteome fraction of tested protein", side = 1, line = 4, outer = TRUE, cex = 1.3)
mtext("Cost / Benefit", side = 2, line = 4, outer = TRUE, cex = 1.3)



dev.off()

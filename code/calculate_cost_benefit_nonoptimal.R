library(here)
library(xtable)

directory <- paste0(here(), "/code")
setwd(directory)

modelname="A7simple"
is.reversible <- 0
predict.parameters <- 0

source("uni_colors.R")
source("initialize_model.R")

opt_data <- read.csv(paste0("../data/", modelname, is.reversible, "_protein_cost.csv"), row.names = 1)
opt_data <- opt_data[opt_data$convergence == 4, ]

row <- 1
rho_cond <- rho_cond[1]
a_cond <- a_cond[,1] #, drop=FALSE]

rows <- length(unique(opt_data$reaction)) *
  nrow(opt_data) *
  length(grep("^v\\.", colnames(opt_data)))

out <- vector("list", rows)
idx <- 1

v_idx <- grep("^v\\.", colnames(opt_data))
f_idx <- grep("^f\\.", colnames(opt_data))
c_idx <- grep("^c\\.", colnames(opt_data))
proteins <- split(opt_data, opt_data$reaction)


for (tested_protein in names(proteins)) {
  one_prot <- proteins[[tested_protein]]
  prot_by_phi <- split(one_prot, one_prot$phi)
  
  for (phi in names(prot_by_phi)) {
    one_phi <- prot_by_phi[[phi]]
    
    vs <- one_phi[row, v_idx]
    fs <- one_phi[row, f_idx]
    cint <- one_phi[row, c_idx]
    cint <- cint[, (1+length(a_cond)):length(cint)] # internal metabolites only

    fint <- cint/rho
    taus <- tau(a_cond, t(cint))
    dtaus <- dtau(a_cond, t(cint))
    transport_benefit_term <- unlist(vs) %*% dtaus %*% unlist(fint)
    growth_rate <- one_phi[row, "mu"]
    M_colSums <- colSums(M)
    
    for(j in 1:length(vs)){
      Mjp <- M["P",j]
      local_cost <- growth_rate*taus[j]
      local_benefit <- unlist(vs) %*% (dtaus %*% M[, j])
      transport_benefit <- M_colSums[j] * transport_benefit_term

      out[[idx]] <- list(
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
      idx <- idx + 1
    }
  }
}
results <- do.call(rbind, lapply(out, as.data.frame))


xlim <- c(0, 1)
ylim <- c(-2.1,2.1)
plots <- c("sum", "local_cost", "local_benefit", "transport_benefit")
plot_names <- c("sum", "(marginal) protein investment",
                "(marginal) local protein value", "biomass production cost")
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
       pch = 19, type = "l", lwd = 3, #  cex = 0.6,
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
    for(plot_idx in length(plots):1){
      # points(one_prot$phi, one_prot[[plots[plot_idx]]],
      #        pch = 19, col = colors[plot_idx], cex=0.6)
      lines(one_prot$phi, one_prot[[plots[plot_idx]]],
             lwd = 3, col = colors[plot_idx])
    }
    
    # only y-axis on first column
    if(prot_idx == 1){
      axis(2, at = c(-2,-1,0,1,2))
      mtext(gsub("v\\.", "", reaction), side = 2, line = 2.5, cex = 0.9)
    }
    
    # only x-axis on last row
    if(rxn_idx == n_prot){
      axis(1)
      mtext(prot, side = 1, line = 2.5, cex = 0.9)
    }
    
    # add legend once in outer margin
    if(rxn_idx == 1 & prot_idx == 1){
      par(xpd = NA)  # allow legend outside plot region
      legend(0, 3.6, horiz = TRUE,
             legend = rev(plot_names), 
             col = rev(colors), 
             pch = 19, bty = "n", cex = 1.25)
      par(xpd = FALSE)
      
    }
  }
}

# Add one common x and y label
mtext(bquote(Phi * " of under-/overexpressed protein"), side = 1, line = 4, outer = TRUE, cex = 1.3)
mtext("Protein assessed for costs and benefits", side = 2, line = 4, outer = TRUE, cex = 1.3)


dev.off()

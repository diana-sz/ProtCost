library(here)

modelname <- "M8"
is.reversible <- 0
predict.parameters <- 0

source(here("code", "uni_colors.R"))
source(here("code", "initialize_model.R"))

opt_data <- read.csv(here("data", paste0(modelname, "_protein_cost.csv")), row.names = 1)
opt_data <- opt_data[opt_data$convergence == 4, ]

row <- 1
rho_cond <- rho_cond[1]
a_cond <- a_cond[,1, drop=FALSE]

rows <- length(unique(opt_data$reaction)) *
  nrow(opt_data) *
  length(grep("^v\\.", colnames(opt_data)))

out <- vector("list", rows)
idx <- 1

v_idx <- grep("^v\\.", colnames(opt_data))
f_idx <- grep("^f\\.", colnames(opt_data))
c_idx <- grep("^c\\.", colnames(opt_data))
proteins <- split(opt_data, opt_data$reaction)


cex_labels <- 0.85

M_colSums <- colSums(M)
for (tested_protein in unique(opt_data$reaction)) {
  one_prot <- proteins[[tested_protein]]
  prot_by_phi <- split(one_prot, one_prot$phi)
  
  for (phi in names(prot_by_phi)) {
    one_phi <- prot_by_phi[[phi]]
    
    vs <- one_phi[row, v_idx]
    fs <- one_phi[row, f_idx]
    cint <- one_phi[row, c_idx]
    cint <- cint[, (1+length(a_cond)):length(cint)] # internal metabolites only

    fint <- cint/rho_cond
    taus <- tau(a_cond, t(cint))
    dtaus <- dtau(a_cond, t(cint))
    transport_benefit_term <- unlist(vs) %*% dtaus %*% unlist(fint)
    growth_rate <- one_phi[row, "mu"]
    
    for(j in 1:length(vs)){
      Mjp <- M["P",j]
      catalytic_cost <- growth_rate*taus[j]
      saturation_value <- unlist(vs) %*% (dtaus %*% M[, j])
      crowding_value <- M_colSums[j] * transport_benefit_term

      out[[idx]] <- list(
        growth_rate = growth_rate,
        tested_protein = tested_protein,
        phi = phi,
        rel_phi = one_phi$rel_phi,
        variable = colnames(vs)[j],
        protein_benefit = as.numeric(Mjp),
        catalytic_cost = -as.numeric(catalytic_cost),
        saturation_value = -as.numeric(saturation_value),
        crowding_value = as.numeric(crowding_value),
        sum = as.numeric(Mjp - catalytic_cost - saturation_value + crowding_value)
      )
      idx <- idx + 1
    }
  }
}
results <- do.call(rbind, lapply(out, as.data.frame))


xlim <- c(0, 4)
ylim <- c(-2.1,2.1)
plots <- c("sum", "catalytic_cost", "saturation_value", "crowding_value")
plot_names <- c("marginal value", "catalytic cost", "saturation value", "crowding value")
colors <- c("black", uni_red, uni_green, uni_lila)

tested_proteins <- unique(results$tested_protein)
variables <- unique(results$variable)
n_prot <- length(unique(tested_proteins))
n_rxn  <- length(variables)

png(here("figures", paste0(modelname, is.reversible, "_cost_benefit.png")), 
    type="cairo", units="cm",
    width=22, height=22, res=300)

par(mfrow = c(n_rxn + 1, n_prot), mar = c(1.7,0,0,0.5), oma = c(5.2,5.7,0.5,0.5), xpd = FALSE)

# plot growth rates on top
prot_idx <- 0
for(protein in tested_proteins){
  prot_idx <- prot_idx + 1
  one_protein <- results[results$tested_protein == protein, ]
  plot(one_protein$rel_phi, one_protein$growth_rate, 
       axes = FALSE, 
       pch = 19, cex = 0.55, # type = "l", lwd = 3
       ylim = c(0, max(one_protein$growth_rate)*1.05), xlim = xlim)
  abline(v=1, #one_protein$phi[which.max(one_protein$growth_rate)], 
         col="grey70", lty=2)
  if(prot_idx == 1){
    axis(2, at = c(0, 0.5, 1))
    mtext("Growth rate", side = 2, line = 2.5, cex = cex_labels)
  }
  box()
}

# plot cost-benefit curves
par(mar = c(0,0,0.5,0.5))
rxn_idx <- 0
for(reaction in rev(variables)){
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
    abline(v=1,
            # one_prot$phi[which.max(one_prot$growth_rate)], 
           col="grey70", lty=2)
    box()
    
    # add all curves
    for(plot_idx in length(plots):1){
      points(one_prot$rel_phi, one_prot[[plots[plot_idx]]],
             pch = 19, col = colors[plot_idx], cex=0.55)
      # lines(one_prot$phi, one_prot[[plots[plot_idx]]],
      #        lwd = 1, col = colors[plot_idx])
    }
    
    # only y-axis on first column
    if(prot_idx == 1){
      axis(2, at = c(-2,-1,0,1,2))
      mtext(gsub("v\\.", "", reaction), side = 2, line = 2.5, cex = 0.75)
    }
    
    # only x-axis on last row
    if(rxn_idx == n_prot){
      axis(1)
      mtext(prot, side = 1, line = 2.5, cex = 0.75)
    }
    
    # add legend once in outer margin
    if(rxn_idx == 1 & prot_idx == 1){
      par(xpd = NA)  # allow legend outside plot region
      legend(0, 3.6, horiz = TRUE,
             legend = rev(plot_names), 
             col = rev(colors), 
             pch = 19, bty = "n", cex = 1.2)
      par(xpd = FALSE)
      
    }
  }
}

# Add one common x and y label
mtext("Proteome fraction of under-/overexpressed protein relative to optimum",
      side = 1, line = 4.4, outer = TRUE, cex = cex_labels)
mtext("Protein assessed for costs and benefits",
      side = 2, line = 4.3, outer = TRUE, cex = cex_labels)


dev.off()

library(here)

modelname <- "M8"
is.reversible <- 0
predict.parameters <- 0

source(here("code", "uni_colors.R"))
source(here("code", "initialize_model.R"))

opt_data <- read.csv(here("data", paste0(modelname, "_protein_cost.csv")))
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
        marginal_value = as.numeric(Mjp - catalytic_cost - saturation_value + crowding_value)
      )
      idx <- idx + 1
    }
  }
}
results <- do.call(rbind, lapply(out, as.data.frame))

# --- numerical derivative of growth_rate w.r.t. rel_phi, per protein ---
# (growth_rate doesn't depend on `variable`, so compute it once per tested_protein
#  rather than per protein/variable combo, to avoid redundant/duplicate work)
results <- results[order(results$tested_protein, results$variable, results$rel_phi), ]
results$derivative <- NA_real_

for (protein in unique(results$tested_protein)) {
  for (variable in unique(results$variable)) {
    idx <- which(results$tested_protein == protein & results$variable == variable)
    if (length(idx) > 1) {
      x <- results$rel_phi[idx]
      y <- results$growth_rate[idx]
      n <- length(x)
      
      dydx <- numeric(n)
      dydx[1] <- (y[2] - y[1]) / (x[2] - x[1])                     # forward diff at start
      dydx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])             # backward diff at end
      if (n > 2) {
        dydx[2:(n - 1)] <- (y[3:n] - y[1:(n - 2)]) / (x[3:n] - x[1:(n - 2)])  # central diff
      }
      
      results$derivative[idx] <- dydx
    }
  }
}


xlim <- c(0, 4)
ylim <- c(-2.1,2.1)
plots <- c("marginal_value", "derivative")
plot_names <- c("marginal value", "derivative of growth rate")
colors <- c("black", uni_red)

tested_proteins <- unique(results$tested_protein)
n_prot <- length(tested_proteins)

png(here("figures", paste0(modelname, "_marginal_val.png")), 
    type="cairo", units="cm",
    width=22, height=8, res=300)

par(mfrow = c(2, n_prot), mar = c(1.7,0,0,0.5), oma = c(5.2,5.7,0.5,0.5), xpd = FALSE)

# --- row 1: growth rate per protein ---
prot_idx <- 0
for(protein in tested_proteins){
  prot_idx <- prot_idx + 1
  one_protein <- results[results$tested_protein == protein, ]
  plot(one_protein$rel_phi, one_protein$growth_rate, 
       axes = FALSE, 
       pch = 19, cex = 0.55,
       ylim = c(0, max(one_protein$growth_rate)*1.05), xlim = xlim)
  abline(v=1, col="grey70", lty=2)
  if(prot_idx == 1){
    axis(2, at = c(0, 0.5, 1))
    mtext("Growth rate", side = 2, line = 2.5, cex = cex_labels)
  }
  box()
}

# --- row 2: marginal value + derivative, only for each protein's own reaction ---
par(mar = c(0,0,0.5,0.5))
prot_idx <- 0
for(protein in tested_proteins){
  prot_idx <- prot_idx + 1
  
  # keep only the variable that matches this protein's own reaction (strip "v." prefix to compare)
  one_prot <- results[results$tested_protein == protein &
                        gsub("^v\\.", "", results$variable) == protein, ]
  
  # plot background
  plot(NA, xlim = xlim, ylim = ylim,
       axes = FALSE, xlab = "", ylab = "")
  abline(h=0)
  abline(v=1, col="grey70", lty=2)
  box()
  
  # add both curves
  for(plot_idx in length(plots):1){
    points(one_prot$rel_phi, one_prot[[plots[plot_idx]]],
           pch = 19, col = colors[plot_idx], cex=0.55)
  }
  
  # only y-axis on first column
  if(prot_idx == 1){
    axis(2, at = c(-2,-1,0,1,2))
    mtext("Marginal value / derivative", side = 2, line = 2.5, cex = 0.75)
  }
  
  # x-axis on every column since this is the last row
  axis(1)
  mtext(protein, side = 1, line = 2.5, cex = 0.75)
  
  # add legend once in outer margin
  if(prot_idx == 1){
    par(xpd = NA)
    legend(0, 3.3, horiz = TRUE,
           legend = plot_names, 
           col = colors, 
           pch = 19, bty = "n", cex = 1.2)
    par(xpd = FALSE)
  }
}

# Add one common x label (dropped the previous y outer label since it
# referenced multiple reaction rows that no longer exist)
mtext(expression("Proteome fraction of under-/overexpressed protein relative to optimum " * italic("\u03A6") * "/" * italic("\u03A6")^"\u204E"),
      side = 1, line = 4.4, outer = TRUE, cex = cex_labels)

dev.off()
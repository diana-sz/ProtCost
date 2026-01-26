library(here)
library(RColorBrewer)

setwd(here())

cex_all <- 1
phi_xlim <- 0.5

modelnames <- c("M10fuel_efflux", "M10fuel_efflux_rev")

for(n in seq_along(modelnames)){
  modelname <- modelnames[n]
  
  results <- read.csv(paste0("data/", modelname, "_tradeoff.csv"))
  results <- results[results$convergence == 4,]  # plot only simulations that converged
  results <- results[results$fuel_phi <= phi_xlim,]
  results$efflux_phi <- results$p.EFFLUX/results$c.P

  # Set colors based on fuel_phi
  palette_colors <- brewer.pal(length(unique(results$kcat)), "Paired")
  color_map <- setNames(palette_colors, unique(results$kcat))
  point_colors <- color_map[as.character(results$kcat)]
  
  png(paste0("figures/", modelname,"_mu.png"),
      title=paste("Protein cost testing", modelname), 
      type="cairo", units="cm", res = 300,
      width=9, height=7.5)

  par(mfrow=c(1,1), oma = c(0,0,0,0))
  par(mar = c(3.7, 3.7, 0.5, 1))
  
  plot(mu ~ fuel_phi, data=results, 
       ylab = NA,
       xlab = NA,
       cex = cex_all,
       pch = 19,
       col = point_colors, 
       xlim = c(0, phi_xlim),
       ylim = c(0, max(results$mu)),
       axes = FALSE)
  
  box()
  axis(1)
  axis(2, las=2)

  legend("topright", col = rev(palette_colors), 
         legend = rev(unique(results$kcat)), pch=19, cex = cex_all, ncol = 2,
         title = expression("efflux pump " * k["cat"]), bty = "n")

  mtext(bquote("Growth rate" ~ "[" * h^-1 * "]"), side = 2, outer = FALSE, line = 2.5, cex = cex_all)
  mtext("Proteome fraction fuel reaction", side = 1, outer = FALSE, line = 2.3, cex = cex_all)
  
  # mtext(
  #   paste0("(", letters[1], ")"),
  #   side = 3,
  #   adj = 0.5,
  #   line = 0.35,
  #   cex = cex_all,
  #   font = 2
  # )
  
  dev.off()
  

  png(paste0("figures/", modelname,"_pump.png"),
      title=paste("Protein cost testing", modelname), 
      type="cairo", units="cm", res = 300,
      width=9, height=7.5)
  
  par(mfrow=c(1,1), oma = c(0,0,0,0))
  par(mar = c(3.7, 3.7, 0.5, 1))
  
  plot(efflux_phi ~ fuel_phi, data=results,
       ylab = NA,
       xlab = NA,
       ylim = c(0, max(results$efflux_phi)),
       xlim = c(0, phi_xlim),
       cex = cex_all,
       pch = 19,
       axes = FALSE,
       col = point_colors)
  
  box()
  axis(1)
  axis(2, las=2) #, at = seq(0, 0.25, 0.05), labels = c(0.0, NA, 0.1, NA, 0.2, NA))

  mtext("Proteome fraction pump", side = 2, outer = FALSE, line = 2.8, cex = cex_all)
  mtext("Proteome fraction fuel reaction", side = 1, outer = FALSE, line = 2.3, cex = cex_all)
  
  # mtext(
  #   paste0("(", letters[2], ")"),
  #   side = 3,
  #   adj = 0.5,
  #   line = 0.35,
  #   cex = cex_all,
  #   font = 2
  # )
  
  dev.off()
  
}



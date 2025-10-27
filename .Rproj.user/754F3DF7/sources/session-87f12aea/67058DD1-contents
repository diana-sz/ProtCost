library(here)
library(RColorBrewer)

setwd(here())

modelnames <- c("A9fuel_efflux", "A9fuel_efflux_rev")

for(modelname in modelnames){
  
  results <- read.csv(paste0("data/", modelname, "_tradeoff_test3.csv"))
  results <- results[results$convergence == 4,]  # plot only simulations that converged
  results$efflux_phi <- results$p.EFFLUX/results$c.p
  
  
  # Set colors based on fuel_phi
  palette_colors <- brewer.pal(length(unique(results$kcat)), "Dark2")
  color_map <- setNames(palette_colors, unique(results$kcat))
  point_colors <- color_map[as.character(results$kcat)]
  
  png(paste0("figures/", modelname,"_tradeoff_test.png"),
      title=paste("Protein cost testing", modelname), 
      type="cairo", units="cm", res = 300,
      width=18, height=8)
  
  par(mar=c(1,4.5,1,1), mfrow = c(1,2), oma = c(3,0,0,0))
  
  plot(mu ~ fuel_phi, data=results, 
       #xlab = "Proteome fraction fuel reaction",
       ylab = bquote("Growth rate" ~ "[" * h^-1 * "]"),
       cex = 1.,
       pch = 19,
       cex.lab = 1.,
       col = point_colors, 
       xlim = c(0,0.8))
  legend("topright", col = palette_colors, 
         legend = unique(results$kcat), pch=19, cex = 0.9,
         title = expression("efflux pump " * k["cat"]), bty = "n")
  plot(efflux_phi ~ fuel_phi, data=results,
       #xlab = "Proteome fraction fuel reaction",
       ylab = "Proteome fraction efflux pump", #  bquote("Growth rate" ~ "[" * h^-1 * "]"),
       ylim = c(0, max(results$efflux_phi)), xlim = c(0,0.8),
       cex = 1.,
       pch = 19,
       cex.lab = 1.,
       col = point_colors)

  
  mtext("Proteome fraction fuel reaction", side = 1, outer = TRUE, line = 1.3, cex = 1)
  
  
  dev.off()
  
}



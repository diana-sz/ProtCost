library(here)
library(RColorBrewer)

cex_lab <- 0.75
xlim <- c(0, 1)
ylim <- c(0, 1)
title_line <- 0.4

plot_composition <- function(proteome, target_phi,
                             colors,
                             xlim,
                             ylim,
                             main = "Proteome Composition",
                             ylab = "Composition",
                             xlab = "Proteome fraction",
                             legend=TRUE,
                             cex_lab=1) {
  # Dimensions
  m <- ncol(proteome)
  
  # Cumulative sum by row (used for stacking)
  proteome <- proteome/rowSums(proteome)
  y_cum <- t(apply(proteome[ncol(proteome):1], 1, cumsum))
  
  # Base plot (empty)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main,
       cex.lab = cex_lab, axes=FALSE)
  axis(1)
  axis(2, labels=FALSE)
  
  # Bottom line (start from 0)
  y_prev <- rep(0, length(target_phi))
  
  # Draw polygons for each protein group
  for (i in m:1) {
    y_top <- y_cum[, i]
    polygon(
      c(target_phi, rev(target_phi)),
      c(y_prev, rev(y_top)),
      col = colors[i],
      border = NA
    )
  }
  
  box()
  
  # Add legend
  if(legend){
    leg_text <- gsub("c\\.|p\\.", "", colnames(proteome))
    par(xpd=NA)
    legend(0.5, 1.52, legend = leg_text, fill = rev(colors), bty = "n", 
           cex = 0.9, ncol = 5, xjust = 0.5)
    par(xpd=FALSE)
  }
}

rxn <- "IDLE"

for(modelname in c("M9_IDLE", "M9_IDLE_rev")){
  data <- read.csv(here("data", paste0(modelname, ".csv")), row.names = 1)
  data <- data[data$convergence == 4, ]
  data <- data[order(data$phi), ]
  
  png(here("figures", paste0(modelname, ".png")), 
      type="cairo", units="cm",
      width=18, height=5.5, res=300)
  par(mfcol=c(1,3), mar = c(3.7,0.8,3.8,0.3), oma = c(0,1.5,0,0))
  
  proteome <- data[, grep("p\\.", colnames(data))]
  biomass <- data[, grep("c\\.", colnames(data))]
  biomass <- biomass[, -grep("x_", colnames(biomass))]
  mu <- data$mu/max(data$mu)
  plotted_phi <- data$phi
    
  plot(plotted_phi, mu,
       xlim = xlim,
       ylim = c(0, max(mu)),
       type = "l",
       lwd = 3,
       axes = FALSE,
       xlab = NA,
       ylab=NA)
  mtext(bquote("Normalized growth rate"), side = 3, cex = cex_lab, line = title_line)
  axis(2, las = 2)
  axis(1)
  box()

  colors <- brewer.pal(ncol(proteome), "PuBu")
  colors[c(2,3,5)] <- c("#E2F1AF", "#A8577E", "#30011E")
  plot_composition(proteome, plotted_phi, colors, main="", 
                   ylab = NA, xlab = NA,
                   xlim = xlim, ylim = ylim,
                   cex_lab=cex_lab,
                   legend=TRUE)
  mtext("Proteome composition", side = 3, cex = cex_lab, line = title_line)
    
  bio_colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
  plot_composition(biomass, plotted_phi, bio_colors, main="", 
                   ylab = NA, xlab = NA,
                   xlim = xlim, ylim = ylim,
                   cex_lab=cex_lab,
                   legend=TRUE)
  mtext("Biomass composition", side = 3, cex = cex_lab, line = title_line)

  mtext(bquote("Idle protein proteome fraction"), side = 1,
        outer = TRUE, cex = cex_lab, line = -1.1, adj = 0.52)
    
  dev.off()
}

library(here)
library(RColorBrewer)
library(dplyr)

cex_lab <- 0.75
ylim <- c(0, 1.05)
title_line <- 0.4

scott <- read.csv(here("data", "experimental", "scott_4B.csv"))
for(m in unique(scott$medium)){
  idx <- scott$medium == m
  one_cond <- scott[idx,]
  one_cond$mu <- one_cond$mu / mean(one_cond$mu[one_cond$b_gal < 1])  # no values at exactly zero
  scott$mu[idx] <- one_cond$mu
}
scott$b_gal <- scott$b_gal/100

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
  axis(2, las = 2) #, labels=FALSE)
  
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
    x_pos <- ifelse(length(leg_text) > 8, -0.25, -0.05)
    legend(x_pos, 1.45, legend = leg_text, fill = rev(colors), bty = "n", 
           cex = 0.9, ncol = 5, xjust = 0)
    par(xpd=FALSE)
  }
}

rxn <- "IDLE"

for(modelname in c("M9_IDLE_rev", "M10_Q_IDLE_rev")){
  data <- read.csv(here("data", paste0(modelname, ".csv")), row.names = 1)
  data <- data[data$convergence == 4, ]
  data <- data[order(data$phi), ]
  max_x <- ifelse(grepl("M9", modelname), 1, 0.45)
  xlim <- c(0, max_x)
  
  
  png(here("figures", paste0(modelname, ".png")), 
      type="cairo", units="cm",
      width=19, height=5.5, res=300)
  par(mfcol=c(1,3), mar = c(0.8,3.9,0.8,0.5), oma = c(3.2,1,2.5,0))
  
  proteome <- data[, grep("p\\.", colnames(data))]
  
  if("p.Q" %in% colnames(proteome)){
    proteome <- proteome %>%
      relocate("p.Q", .after = last_col())
  }

  biomass <- data[, grep("c\\.", colnames(data))]
  biomass <- biomass[, -grep("x_", colnames(biomass))]
  mu <- data$mu/max(data$mu)
  plotted_phi <- data$phi
    
  plot(plotted_phi, mu,
       xlim = xlim,
       ylim = ylim,
       type = "l",
       lwd = 3,
       axes = FALSE,
       xlab = NA,
       ylab=NA)
  
  if("p.Q" %in% colnames(proteome)){
    points(scott$b_gal, scott$mu, col = "grey50", pch = 24)
  }
  #mtext(bquote("Normalized growth rate"), side = 3, cex = cex_lab, line = title_line)
  mtext(bquote("Growth rate relative to optimum " * mu / mu^"\u204E"), side = 2, cex = cex_lab*0.9, line = 2.3)
  
  axis(2, las = 2)
  axis(1)
  box()

  colors <- brewer.pal(9, "PuBu")
  colors[c(2,3,5)] <- c("#E2F1AF", "#A8577E", "#30011E")
  
  if("p.Q" %in% colnames(proteome)){
    colors <- c("grey80", colors)}
  plot_composition(proteome, plotted_phi, colors, main="", 
                   ylab = NA, xlab = NA,
                   xlim = xlim, ylim = ylim,
                   cex_lab = cex_lab,
                   legend = TRUE)
  #mtext("Proteome composition", side = 3, cex = cex_lab, line = title_line)
  mtext("Proteome mass fraction", side = 2, cex = cex_lab*0.9, line = 2.3)
  
    
  bio_colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
  plot_composition(biomass, plotted_phi, bio_colors, main="", 
                   ylab = NA, xlab = NA,
                   xlim = xlim, ylim = ylim,
                   cex_lab=cex_lab,
                   legend=TRUE)
  #mtext("Biomass composition", side = 3, cex = cex_lab, line = title_line)
  mtext("Biomass fraction", side = 2, cex = cex_lab*0.9, line = 2.3)
  
  mtext(bquote("Proteome fraction of idle protein " * italic("\u03A6")["IDLE"]), side = 1,
        outer = TRUE, cex = cex_lab, line = 1.8, adj = 0.52)
    
  dev.off()
}

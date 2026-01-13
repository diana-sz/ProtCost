library(here)
library(RColorBrewer)

setwd(here())

cex_lab <- 0.75
xlim <- c(0, 1)
ylim <- c(0, 1)
xlab_line <- 1.6
title_line <- 0.4

plot_composition <- function(proteome, target_phi, colors = NULL,
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
  plot(NA, xlim = xlim, ylim =ylim, xlab = xlab, ylab = ylab, main = main,
       cex.lab = cex_lab)
  
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
  
  # Add legend
  if(legend){
    leg_text <- gsub("c\\.", "", colnames(proteome))
    leg_text <- gsub("p\\.", "", leg_text)
    par(xpd=NA)
    legend(0.45, 1.45, legend = leg_text, fill = rev(colors), bty = "n",
           cex = 0.82, ncol = ceiling(m/2), xjust = 0.5)
    par(xpd=FALSE)
  }
}

rxn <- "GFP"

for(modelname in c("M9_GFP", "M9_GFP_rev")){
  data <- read.csv(paste0("data/", modelname, ".csv"), row.names = 1)
  data <- data[data$convergence == 4, ]
  data <- data[complete.cases(data),]
  
  png(paste0("figures/", modelname, ".png"), 
      type="cairo", units="cm",
      width=18, height=6, res=300)
  par(mfcol=c(1,3), mar = c(3.7,2,3.8,0.5))
  
  data <- data[order(data$phi), ]
  proteome <- data[, grep("p\\.", colnames(data))]
  biomass <- data[, grep("c\\.", colnames(data))]
  biomass <- biomass[, -grep("x_", colnames(biomass))]
  mu <- data$mu
  plotted_phi <- data$phi
    
  plot(plotted_phi, mu,
       xlim = xlim, ylim = c(0, max(mu)),
       pch=20,
       xlab = NA,
       ylab=NA)
  mtext(bquote("Growth rate"), side = 3, cex = cex_lab, line = title_line)
  
  colors <- brewer.pal(ncol(proteome), "PuBu") #PuBU #YlGnBu
  plot_composition(proteome, plotted_phi, colors, main="", 
                   ylab=NA, #"Proteome composition",
                   xlab = NA,
                   cex_lab=cex_lab,
                   legend=TRUE)
  mtext("Proteome composition", side = 3, cex = cex_lab, line = title_line)
    
  colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
  plot_composition(biomass, plotted_phi, colors, main="", 
                   ylab = NA, # "Biomass composition",
                   xlab = NA,
                   cex_lab=cex_lab,
                   legend=TRUE)
  mtext("Biomass composition", side = 3, cex = cex_lab, line = title_line)
  mtext(bquote("Idle protein proteome fraction"), side = 1,
        outer = TRUE, cex = cex_lab, line = -1.1, adj = 0.52)
    
    
  dev.off()
}





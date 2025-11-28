library(here)
library(RColorBrewer)

setwd(here())
cex_lab <- 0.75
xlim <- c(0, 0.8)
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
  n <- nrow(proteome)
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
    legend(0.4, 1.53, legend = leg_text, fill = rev(colors), bty = "n", cex = 0.82, ncol = 4, xjust = 0.5)
    par(xpd=FALSE)
  }
}



for(modelname in c("A7simple0", "A7simple_rev1")){
  data <- read.csv(paste0("data/", modelname, "_protein_cost.csv"), row.names = 1)
  data <- data[data$convergence == 4, ]

  for(rxn in unique(data$reaction)){
    
    png(paste0("figures/", modelname, "_", rxn, ".png"), 
        type="cairo", units="cm",
        width=18, height=5.2, res=300)
    par(mfcol=c(1,4), mar = c(3.7,2,3.8,0.5))
    
    one_prot <- data[data$reaction == rxn, ]
    one_prot <- one_prot[order(one_prot$phi), ]
    target_phi <- one_prot$phi
    proteome <- one_prot[, grep("p\\.", colnames(one_prot))]
    biomass <- one_prot[, grep("c\\.", colnames(one_prot))]
    biomass <- biomass[, -grep("x_", colnames(biomass))]
    mu <- one_prot$mu
    
    plot(one_prot$phi/one_prot$phi[which.max(one_prot$mu)], mu,
         xlim = c(0,  4), ylim = c(0, max(mu)),
         pch=20,
         xlab = NA,
         ylab=NA)
    text(3.2, 0.2, rxn, cex=cex_lab*1.5)
    mtext(bquote("Growth rate"), side = 3, cex = cex_lab, line = title_line)
    mtext(bquote(Phi * "/" * Phi[opt]), side = 1, 
          cex = cex_lab*1.1, line = xlab_line, padj = 1)
    
    plot(one_prot$phi, mu, xlim = xlim, ylim = c(0, max(mu)),
         pch=20,
         #main = rxn,
         #cex.lab = cex_lab*1.3,
         xlab = NA, # "Proteome fraction", 
         ylab=NA) #bquote("Growth rate" ~ "[" * h^-1 * "]"))
    mtext(bquote("Growth rate"), side = 3, cex = cex_lab, line = title_line)
    mtext(bquote(Phi), side = 1, 
          cex = cex_lab*1.1, line = xlab_line, padj = 1)
    
    colors <- brewer.pal(ncol(proteome), "PuBu") #PuBU #YlGnBu
    plot_legend <- ifelse(rxn == "TC", TRUE, FALSE)
    plot_composition(proteome, target_phi, colors, main="", 
                     ylab=NA, #"Proteome composition",
                     xlab = NA,
                     cex_lab=cex_lab,
                     legend=plot_legend)
    mtext("Proteome composition", side = 3, cex = cex_lab, line = title_line)
    mtext(bquote(Phi), side = 1, 
          cex = cex_lab*1.1, line = xlab_line, padj = 1)
    
    colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
    plot_composition(biomass, target_phi, colors, main="", 
                     ylab = NA, # "Biomass composition",
                     xlab = NA,
                     cex_lab=cex_lab,
                     legend=plot_legend)
    mtext("Biomass composition", side = 3, cex = cex_lab, line = title_line)
    mtext(bquote(Phi), side = 1, 
          cex = cex_lab*1.1, line = xlab_line, padj = 1)
    
    # 
    # mtext(bquote("Proteome fraction of tested protein" ~ Phi), side = 1, 
    #       outer = TRUE, cex = cex_lab*1.1, line = -1.1, adj = 0.7)
    
    
    dev.off()
  }
}




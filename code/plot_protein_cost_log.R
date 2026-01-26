library(here)
library(RColorBrewer)

setwd(here())

relative_phi <- FALSE

cex_lab <- 0.75
xlim <- c(0.01, ifelse(relative_phi, 5, 1))
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
  plot(NA, 
       xlim = xlim, 
       ylim =ylim, 
       xlab = xlab, 
       ylab = ylab, 
       main = main,
       log = "x",
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
    legend(0.02, 1.52, legend = leg_text, fill = rev(colors), bty = "n", 
           cex = 0.9, ncol = 4, xjust = 0)
    par(xpd=FALSE)
  }
}


for(modelname in c("B", "B_rev")){#  c("M8b0", "M8b_rev1")){ #A7simple_rev1
  data <- read.csv(paste0("data/", modelname, "_protein_cost.csv"), row.names = 1)
  data <- data[data$convergence == 4, ]
  data <- data[complete.cases(data),]
  data <- data[data$phi > 0.01, ]
  
  for(rxn in unique(data$reaction)){
    suffix <- ifelse(relative_phi, "_rel", "_abs")
    png(paste0("figures/", modelname, "_", rxn, suffix, "_log.png"), 
        type="cairo", units="cm",
        width=20, height=5.5, res=300)
    par(mfcol=c(1,3), mar = c(3.7,2,3.8,0.5))
    
    one_prot <- data[data$reaction == rxn, ]
    one_prot <- one_prot[order(one_prot$phi), ]
    proteome <- one_prot[, grep("p\\.", colnames(one_prot))]
    q_col <- which(colnames(proteome)=="p.Q")
    if(length(q_col) == 1){
      proteome <- cbind(proteome[, q_col], proteome[, -q_col])
      colnames(proteome)[1] <- "Q"
    }
    
    biomass <- one_prot[, grep("c\\.", colnames(one_prot))]
    biomass <- biomass[, -grep("x_", colnames(biomass))]
    mu <- one_prot$mu
    one_prot$rel_phi <- one_prot$phi/one_prot$phi[which.max(one_prot$mu)]
    
    plotted_phi <- one_prot$phi
    type <- "p"
    if(relative_phi){
      plotted_phi <- one_prot$rel_phi
      type <- "l"
    }
    
    plot(plotted_phi, mu,
         xlim = xlim,
         ylim = c(0, max(mu)),
         type = type, lwd = 3,
         log = "x",
         cex = 0.6,
         pch=20,
         xlab = NA,
         ylab=NA)
    text(ifelse(relative_phi, 3.2, 0.6), 0.2, rxn, cex=cex_lab*1.5)
    mtext(bquote("Growth rate"), side = 3, cex = cex_lab, line = title_line)
    abline(v = plotted_phi[which.max(mu)], lty = 2, col = "grey70")
    
    colors <- brewer.pal(ncol(proteome), "PuBu") #PuBU #YlGnBu
    plot_legend <- ifelse(rxn == "TS", TRUE, FALSE)
    plot_composition(proteome, plotted_phi, colors, main="", 
                     ylab=NA, #"Proteome composition",
                     xlab = NA,
                     cex_lab=cex_lab,
                     legend=plot_legend)
    mtext("Proteome composition", side = 3, cex = cex_lab, line = title_line)

    colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
    plot_composition(biomass, plotted_phi, colors, main="", 
                     ylab = NA,
                     xlab = NA,
                     cex_lab=cex_lab,
                     legend=plot_legend)
    mtext("Biomass composition", side = 3, cex = cex_lab, line = title_line)

    xlabel <- "Proteome fraction"
    if(relative_phi){
      xlabel <- "Proteome fraction relative to optimum"
    }
    
    mtext(xlabel, side = 1,
          outer = TRUE, cex = cex_lab, line = -1.1, adj = 0.52)
    
    
    dev.off()
  }
}




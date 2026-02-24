library(here)
library(RColorBrewer)

relative_phi <- TRUE
models_to_plot <- c("M8", "M8_rev", "M9_Q", "M9_Q_rev", "B", "B_rev")

cex_lab <- 0.75
xlim <- c(0, ifelse(relative_phi, 5, 1))
ylim <- c(0, 1)
title_line <- 0.4

reaction_names <- list("TS" = "Transporter",
                       "ADPS" = "ADP synthesis",
                       "ATPS" = "ATP synthesis", 
                       "AAS" = "Amino acid synthesis",
                       "NTS" = "Nucleotide synthesis",
                       "RNAP" = "RNA polymerase",
                       "DNAP" = "DNA polymerase",
                       "R" = "Ribosome")


plot_composition <- function(proteome, target_phi,
                             colors,
                             xlim, ylim,
                             main = NA,
                             ylab = NA, xlab = NA,
                             legend = TRUE,
                             cex_lab = 1) {
  
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
  for (i in ncol(proteome):1) {
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
    n <- length(leg_text)
    ncol_legend <- ifelse(n <= 3, n, ceiling(n / 2))
    legend(0, 1.52, legend = leg_text, fill = rev(colors), bty = "n", 
           cex = 0.9, ncol = ncol_legend, xjust = 0)
    par(xpd=FALSE)
  }
}


for(modelname in models_to_plot){
  data <- read.csv(here("data", paste0(modelname, "_protein_cost.csv")), row.names = 1)
  data <- data[data$convergence == 4, ]
  print(paste(modelname, max(data$mu)))

  for(rxn in unique(data$reaction)){
    suffix <- ifelse(relative_phi, "_rel", "_abs")
    png(here("figures", paste0(modelname, "_", rxn, suffix, ".png")),
        type="cairo", units="cm",
        width=18, height=5.5, res=300)
    par(mfcol=c(1,3), mar = c(3.7,0.8,3.8,0.3), oma = c(0,1.5,0,0))
    
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
    mu <- one_prot$mu_norm
    one_prot$rel_phi <- one_prot$phi/one_prot$phi[which.max(one_prot$mu)]
    
    plotted_phi <- one_prot$phi
    plotted_phi_mu <- plotted_phi
    type <- "p"
    if(relative_phi){
      plotted_phi_mu <- c(0, one_prot$rel_phi) # 0 added for visualisation
      plotted_phi <- one_prot$rel_phi
      mu <- c(0, mu)
      type <- "l"
    }
    
    plot(plotted_phi_mu, mu,
         xlim = xlim, ylim = c(0, max(mu)),
         type = type,
         lwd = 3,
         cex = 0.6,
         pch = 20,
         axes = FALSE,
         xlab = NA, ylab = NA)
    text(ifelse(relative_phi, 3.2, 0.7), 0.25, rxn, cex=cex_lab*1.4)
    text(ifelse(relative_phi, 3.2, 0.7), 0.1, reaction_names[[rxn]], cex=cex_lab*1.2)
    axis(2, las = 2)
    axis(1)
    box()
    
    mtext(bquote("Normalized growth rate"), side = 3, cex = cex_lab, line = title_line)
    abline(v = plotted_phi_mu[which.max(mu)], lty = 2, col = "grey70")
    
    colors <- brewer.pal(ncol(proteome), "PuBu")
    if(grepl("M", modelname)){
      colors[c(2,4)] <- c( "#A8577E", "#30011E") # add some contrasting colors
    }
    plot_legend <- ifelse(rxn == "TS", TRUE, FALSE)
    plot_composition(proteome, plotted_phi, colors, main="", 
                     ylab = NA, xlab = NA,
                     xlim = xlim, ylim =ylim,
                     cex_lab = cex_lab,
                     legend = plot_legend)
    mtext("Proteome composition", side = 3, cex = cex_lab, line = title_line)

    colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
    plot_composition(biomass, plotted_phi, colors, main="", 
                     ylab = NA, xlab = NA,
                     xlim = xlim, ylim =ylim,
                     cex_lab = cex_lab,
                     legend = plot_legend)
    mtext("Biomass composition", side = 3, cex = cex_lab, line = title_line)

    xlabel <- "Proteome fraction"
    if(relative_phi){
      xlabel <- "Proteome fraction relative to optimum"
    }
    
    mtext(xlabel, side = 1,
          outer = TRUE, cex = cex_lab, line = -1.2, adj = 0.51)
    
    dev.off()
  }
}


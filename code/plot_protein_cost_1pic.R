library(here)
library(RColorBrewer)

relative_phi <- TRUE
models_to_plot <- c("M8", "M8_rev", "M9_Q")#, "B", "B_rev")

cex_lab <- 0.75
xlim <- c(0, ifelse(relative_phi, 5, 1))
ylim <- c(0, 1)
title_line <- 0.4
axis_label_line <- 2.4

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
  axis(1, labels = axis_labels)
  axis(2, las = 2)
  
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
    legend(-0.03*max(target_phi), 1.48, legend = leg_text, fill = rev(colors), bty = "n", 
           cex = 0.95, ncol = ncol_legend, xjust = 0)
    par(xpd=FALSE)
  }
}


for(modelname in models_to_plot){
  data <- read.csv(here("data", paste0(modelname, "_protein_cost.csv")))
  data <- data[data$convergence == 4, ]
  
  if(modelname == "M8_rev"){
    data_irrev <- read.csv(here("data", "M8_protein_cost.csv"))
    data_irrev <- data_irrev[data_irrev$convergence == 4, ]
  }
  
  print(paste(modelname, max(data$mu)))
  n_rxns <- 8 #length(unique(data$reaction))
  
  suffix <- ifelse(relative_phi, "_rel", "_abs")
  png(here("figures", paste0(modelname, suffix, ".png")),
      type="cairo", units="cm",
      width=22.5, height=3.5*n_rxns, res=300)
  par(mfrow=c(n_rxns, 3), mar = c(0.8,4.5,0.8,0.5), oma = c(3.3,1,3,0))
  
  for(rxn in unique(data$reaction)){
    if(rxn == "Q"){next}
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
    mu <- one_prot$mu/max(one_prot$mu)
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
    
    if(modelname == "M8_rev"){
      one_prot_irrev <- data_irrev[data_irrev$reaction == rxn, ]
      one_prot_irrev <- one_prot_irrev[order(one_prot_irrev$phi), ]
      mu_irrev <- one_prot_irrev$mu/max(one_prot_irrev$mu)
      plotted_phi_mu_irr <- one_prot_irrev$phi/one_prot_irrev$phi[which.max(one_prot_irrev$mu)]
      if(relative_phi & (rxn %in% c("TS", "ATPS", "ADPS", "NTS", "AAS"))){
        plotted_phi_mu_irr <- c(0, plotted_phi_mu_irr) # 0 added for visualisation
        mu_irrev <- c(0, mu_irrev)
        lines(plotted_phi_mu_irr, mu_irrev, col = "grey70")
      }

      
    }
    text(ifelse(relative_phi, 2.65, 0.7), 0.25, rxn, cex=cex_lab*1.4)
    text(ifelse(relative_phi, 2.65, 0.7), 0.1, reaction_names[[rxn]], cex=cex_lab*1.2)
    axis(2, las = 2)
    
    axis_labels <- rxn == "R"
    axis(1, labels = axis_labels)
  
    box()
    
    if(rxn == "NTS"){
      #mtext(bquote("Normalized growth rate"), side = 3, cex = cex_lab, line = title_line)
    # }
    # if(rxn == "NTS"){
      mtext(bquote("Growth rate relative to optimum " * mu / mu^"\u204E"), 
            side = 2, cex = cex_lab, line = axis_label_line, at = -0.2)
    }
    
    opt_phi <- plotted_phi_mu[which.max(mu)]
    abline(v = opt_phi, lty = 2, col = "grey70")
    if(rxn == "ADPS"){print(opt_phi)}
    
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
    
    if(rxn == "NTS"){
      #mtext("Proteome composition", side = 3, cex = cex_lab, line = title_line)
    # }
    # if(rxn == "NTS"){
      mtext("Proteome mass fraction", side = 2, cex = cex_lab, line = axis_label_line, at = -0.2)
    }

    colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
    plot_composition(biomass, plotted_phi, colors, main="", 
                     ylab = NA, xlab = NA,
                     xlim = xlim, ylim =ylim,
                     cex_lab = cex_lab,
                     legend = plot_legend)
    
    if(rxn == "NTS"){
      #mtext("Biomass composition", side = 3, cex = cex_lab, line = title_line)
    # }
    # if(rxn == "NTS"){
      mtext("Biomass fraction", side = 2, cex = cex_lab, line = axis_label_line, at = -0.2)
    }

    xlabel <- expression("Proteome fraction " * italic("\u03A6"))
    if(relative_phi){
      xlabel <- expression("Proteome fraction relative to optimum " * italic("\u03A6") * "/" * italic("\u03A6")^"\u204E")
    }
    mtext(xlabel, side = 1,
          outer = TRUE, cex = cex_lab, line = 1.65, adj = 0.52)
  }
  dev.off()
}


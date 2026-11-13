library(here)
library(RColorBrewer)
source(here("code", "uni_colors.R"))

# Save each figure as TIF (journal submission) and PNG (Overleaf)
open_figure <- function(file_name, width, height, format) {
  if (format == "tif") {
    tiff(here("figures", paste0(file_name, ".tif")),
         type = "cairo", units = "cm", width = width, height = height,
         res = 1000, compression = "lzw", family = "Helvetica")
  } else {
    png(here("figures", paste0(file_name, ".png")),
        type = "cairo", units = "cm", width = width, height = height,
        res = 300, family = "Helvetica")
  }
  par(family = "Helvetica")
}

relative_phi <- TRUE
models_to_plot <- c("M8", "M8_rev", "M9_Q_rev")

cex_lab <- 0.75
xlim <- c(0, ifelse(relative_phi, 5, 1))
ylim <- c(0, 1)
title_line <- 0.4
axis_label_line <- 2.4
# Names of the figure files in the manuscript; other outputs keep their model name
figure_names <- c(M8_rev_rel = "Figure_2",
                  M8_rel = "Appendix_Figure_S7",
                  M8_rev_abs = "Appendix_Figure_S8",
                  M8_abs = "Appendix_Figure_S9",
                  M9_Q_rev_rel = "Appendix_Figure_S10")

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
                             cex_lab = 1,
                             shade_gap = FALSE,
                             shade_from = NULL) {
  
  # Cumulative sum by row (used for stacking)
  proteome <- proteome/rowSums(proteome)
  y_cum <- t(apply(proteome[ncol(proteome):1], 1, cumsum))
  
  # Base plot (empty)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main,
       cex.lab = cex_lab, axes=FALSE, xaxs = "i", yaxs = "i")
  axis(1, labels = axis_labels)
  axis(2, las = 2)
  
  # Shade region beyond the highest observed phi (only when the run
  # stopped converging rather than actually reaching near-zero growth rate)
  if(shade_gap && !is.null(shade_from) && shade_from < xlim[2]){
    rect(shade_from, par("usr")[3], xlim[2], par("usr")[4],
         col = "grey95", border = NA)
  }
  
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
    legend(-0.03*max(target_phi), 1.45, legend = leg_text, fill = rev(colors), bty = "n", 
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
  file_name <- paste0(modelname, suffix)
  if (file_name %in% names(figure_names)) file_name <- figure_names[[file_name]]
  for (format in c("tif", "png")) {
    open_figure(file_name, 22.5, 3.5*n_rxns, format)
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
      }
      
      plot(NA,
           xaxs = "i", yaxs = "i",
           xlim = xlim, ylim = c(0, 1.03),
           axes = FALSE,
           xlab = NA, ylab = NA)
      
      # Shade region beyond the highest observed phi, but only if the run
      # stopped converging rather than actually reaching near-zero growth rate.
      # This decision is reused below for the composition panels too.
      max_phi <- max(plotted_phi_mu, na.rm = TRUE)
      mu_at_max_phi <- mu[which.max(plotted_phi_mu)]
      growth_near_zero_threshold <- 0.1  # mu is normalized so max(mu) == 1
      shade_gap <- max_phi < xlim[2] && mu_at_max_phi > growth_near_zero_threshold
      if(shade_gap){
        rect(max_phi, par("usr")[3], xlim[2], par("usr")[4],
             col = "grey95", border = NA)
      }
      
      
      if(modelname == "M8_rev"){
        one_prot_irrev <- data_irrev[data_irrev$reaction == rxn, ]
        one_prot_irrev <- one_prot_irrev[order(one_prot_irrev$phi), ]
        mu_irrev <- one_prot_irrev$mu/max(one_prot_irrev$mu)
        plotted_phi_mu_irr <- one_prot_irrev$phi/one_prot_irrev$phi[which.max(one_prot_irrev$mu)]
        if(relative_phi & (rxn %in% c("TS", "ATPS", "ADPS", "NTS", "AAS"))){
          plotted_phi_mu_irr <- c(0, plotted_phi_mu_irr) # 0 added for visualisation
          mu_irrev <- c(0, mu_irrev)
          lines(plotted_phi_mu_irr, mu_irrev, col = "grey70", lwd=2)
        }
      }
      
      points(plotted_phi_mu, mu,
             type = "l",
             col = uni_wine,
             lwd = 4)
      
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
      
      colors <- brewer.pal(ncol(proteome), "PuBu")
      if(grepl("M", modelname)){
        colors[c(2,4)] <- c( "#A8577E", "#30011E") # add some contrasting colors
      }
      plot_legend <- ifelse(rxn == "TS", TRUE, FALSE)
      plot_composition(proteome, plotted_phi, colors, main="", 
                       ylab = NA, xlab = NA,
                       xlim = xlim, ylim =ylim,
                       cex_lab = cex_lab,
                       legend = plot_legend,
                       shade_gap = shade_gap,
                       shade_from = max_phi)
      
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
                       legend = plot_legend,
                       shade_gap = shade_gap,
                       shade_from = max_phi)
      
      if(rxn == "NTS"){
        mtext("Biomass fraction", side = 2, cex = cex_lab, line = axis_label_line, at = -0.2)
      }
      
      xlabel <- expression("Proteome fraction " * italic("\u03D5"))
      if(relative_phi){
        xlabel <- expression("Proteome fraction relative to optimum " * italic("\u03D5") * "/" * italic("\u03D5")^"\u204E")
      }
      mtext(xlabel, side = 1,
            outer = TRUE, cex = cex_lab, line = 1.65, adj = 0.52)
    }
    dev.off()
  }
  
  # ---- Source data (Figure 2 only) ----
  if (file_name == "Figure_2") {
    source_dir <- here("figures", "source_data_figures", "Figure 2")
    out_dir <- file.path(source_dir, "2A")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # wine lines and composition panels: reversible kinetics
    sim <- do.call(rbind, lapply(setdiff(unique(data$reaction), "Q"), function(rxn) {
      one <- data[data$reaction == rxn, ]
      one <- one[order(one$phi), ]
      proteome <- one[, grep("p\\.", colnames(one))]
      proteome <- proteome / rowSums(proteome)
      colnames(proteome) <- paste0("proteome_", sub("^p\\.", "", colnames(proteome)))
      biomass <- one[, grep("c\\.", colnames(one))]
      biomass <- biomass[, -grep("x_", colnames(biomass))]
      biomass <- biomass / rowSums(biomass)
      colnames(biomass) <- paste0("biomass_", sub("^c\\.", "", colnames(biomass)))
      data.frame(reaction = rxn,
                 phi = one$phi,
                 phi_relative = one$phi / one$phi[which.max(one$mu)],
                 growth_rate = one$mu,
                 growth_rate_relative = one$mu / max(one$mu),
                 proteome, biomass, check.names = FALSE)
    }))
    sim <- sim[sim$phi_relative <= xlim[2], ]
    write.csv(sim, file.path(out_dir, "Figure_2A_reversible.csv"), row.names = FALSE)
    
    # grey lines: irreversible kinetics
    irr_rxns <- intersect(c("TS", "ATPS", "ADPS", "NTS", "AAS"), unique(data_irrev$reaction))
    irr <- do.call(rbind, lapply(irr_rxns, function(rxn) {
      one <- data_irrev[data_irrev$reaction == rxn, ]
      one <- one[order(one$phi), ]
      data.frame(reaction = rxn,
                 phi = one$phi,
                 phi_relative = one$phi / one$phi[which.max(one$mu)],
                 growth_rate = one$mu,
                 growth_rate_relative = one$mu / max(one$mu))
    }))
    irr <- irr[irr$phi_relative <= xlim[2], ]
    write.csv(irr, file.path(out_dir, "Figure_2A_irreversible.csv"), row.names = FALSE)
    
    writeLines(c(
      "Source data for Figure 2",
      "",
      "Simulations: GBA model M8. Each row of the figure varies the proteome fraction of one reaction.",
      "Only converged simulations with phi_relative <= 5 are included, as in the figure.",
      "The point (0, 0) added to the growth rate lines for display is not included.",
      "",
      "2A/Figure_2A_reversible.csv: reversible kinetics (wine lines, proteome and biomass panels)",
      "  reaction: reaction whose proteome fraction is varied",
      "  phi: proteome fraction of that reaction",
      "  phi_relative: phi relative to its value at the maximal growth rate (phi/phi*)",
      "  growth_rate: growth rate [1/h]",
      "  growth_rate_relative: growth rate relative to its maximum (mu/mu*)",
      "  proteome_<X>: mass fraction of protein X in the proteome (each row sums to 1)",
      "  biomass_<X>: mass fraction of component X in the biomass (each row sums to 1)",
      "",
      "2A/Figure_2A_irreversible.csv: irreversible kinetics (grey lines)",
      "  Columns reaction to growth_rate_relative as above, for TS, ATPS, ADPS, NTS and AAS."
    ), file.path(source_dir, "README.txt"))
  }
}
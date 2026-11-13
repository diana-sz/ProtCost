library(here)
library(RColorBrewer)
library(dplyr)
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
        res = 600, family = "Helvetica")
  }
  par(family = "Helvetica")
}

# Bold panel letter in the top-left corner of the current panel.
# A panel at the edge of the figure also includes the outer margin there.
add_panel_label <- function(label, cex = 1.2) {
  left <- grconvertX(0, "nfc", "ndc")
  top <- grconvertY(1, "nfc", "ndc")
  if (abs(left - grconvertX(0, "nic", "ndc")) < 1e-6) left <- 0
  if (abs(top - grconvertY(1, "nic", "ndc")) < 1e-6) top <- 1
  offset <- 0.02  # inches from the corner
  x <- left + grconvertX(offset, "inches", "ndc")
  y <- top - grconvertY(offset, "inches", "ndc")
  text(grconvertX(x, "ndc", "user"), grconvertY(y, "ndc", "user"),
       label, adj = c(0, 1), font = 2, cex = cex, xpd = NA)
}

cex_lab <- 0.75
ylim <- c(0, 1)
title_line <- 0.4
leg_size <- 0.85

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
       xaxs = "i", yaxs = "i",
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
    x_pos <- ifelse(length(leg_text) > 8, -0.12*ylim[2], -0.02)
    legend(x_pos+0.02, 1.3, legend = leg_text, fill = rev(colors), bty = "n", 
           cex = leg_size, ncol = 5, xjust = 0)
    par(xpd=FALSE)
  }
}

# Write one source data file into figures/source_data_figures/<figure>/<panel>
write_source <- function(df, figure, panel, file_name) {
  out_dir <- here("figures", "source_data_figures", figure, panel)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(out_dir, file_name), row.names = FALSE)
}

rxn <- "IDLE"
# Names of the figure files in the manuscript; other outputs keep their model name
figure_names <- c(M10_Q_IDLE_rev = "Figure_3",
                  M9_IDLE_rev = "Appendix_Figure_S12")

for(modelname in c("M10_Q_IDLE_rev", "M9_IDLE_rev")){
  data <- read.csv(here("data", paste0(modelname, ".csv")), row.names = 1)
  data <- data[data$convergence == 4, ]
  data <- data[order(data$phi), ]
  max_x <- ifelse(grepl("M9", modelname), 1, max(data$phi))
  xlim <- c(0, max_x)
  
  
  file_name <- modelname
  if (file_name %in% names(figure_names)) file_name <- figure_names[[file_name]]
  for (format in c("tif", "png")) {
    open_figure(file_name, 19, 5.5, format)
    par(mfcol=c(1,3), mar = c(0.7,3.9,0.9,0.5), oma = c(3.2,1,2.5,0))
    
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
         xaxs = "i", yaxs = "i",
         col = uni_wine,
         xlim = xlim,
         ylim = c(0, 1.08),
         type = "l",
         lwd = 3,
         axes = FALSE,
         xlab = NA,
         ylab=NA)
    
    if("p.Q" %in% colnames(proteome)){
      points(scott$b_gal, scott$mu, col = "grey50", pch = 24)
    }
    mtext(bquote("Growth rate relative to optimum " * mu / mu^"\u204E"), side = 2, cex = cex_lab*0.89, line = 2.3)
    
    add_panel_label(LETTERS[1], cex = cex_lab*1.6)
    
    if("p.Q" %in% colnames(proteome)){
      par(xpd=NA)
      legend(0, 1.38,
             legend = c("Model", "Experiment (Scott et al. 2010)"),
             lty = c(1, NA),
             lwd = c(3, 1),
             pch = c(NA, 2),
             cex = leg_size,
             col = c(uni_wine, "grey50"),
             bty = "n")
      par(xpd=FALSE)
    }
    
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
    mtext("Proteome mass fraction", side = 2, cex = cex_lab*0.9, line = 2.3)
    
    add_panel_label(LETTERS[2], cex = cex_lab*1.6)
    
    bio_colors <- rev(brewer.pal(ncol(biomass), "RdBu"))
    plot_composition(biomass, plotted_phi, bio_colors, main="", 
                     ylab = NA, xlab = NA,
                     xlim = xlim, ylim = ylim,
                     cex_lab=cex_lab,
                     legend=TRUE)
    mtext("Biomass fraction", side = 2, cex = cex_lab*0.9, line = 2.3)
    
    mtext(bquote("Proteome fraction of idle protein " * italic("\u03D5")["IDLE"]), side = 1,
          outer = TRUE, cex = cex_lab, line = 1.8, adj = 0.52)
    
    add_panel_label(LETTERS[3], cex = cex_lab*1.6)
    
    dev.off()
  }
  
  # ---- Source data (main text figure only) ----
  if (file_name == "Figure_3") {
    fig_dir <- "Figure 3"
    # Only the points inside the plotted axis range
    in_x <- plotted_phi >= xlim[1] & plotted_phi <= xlim[2]
    
    # Panel A: model curve and experimental points
    panel_a_model <- data.frame(phi_IDLE = plotted_phi[in_x],
                                relative_growth_rate = mu[in_x])
    write_source(panel_a_model, fig_dir, "3A", "Figure_3A_model.csv")
    
    scott_in <- scott$b_gal >= xlim[1] & scott$b_gal <= xlim[2] &
      scott$mu >= ylim[1] & scott$mu <= ylim[2]
    panel_a_exp <- data.frame(phi_IDLE = scott$b_gal[scott_in],
                              relative_growth_rate = scott$mu[scott_in])
    write_source(panel_a_exp, fig_dir, "3A", "Figure_3A_Scott_2010.csv")
    
    # Panel B: proteome composition (mass fractions, as plotted)
    proteome_frac <- proteome / rowSums(proteome)
    colnames(proteome_frac) <- sub("^p\\.", "", colnames(proteome_frac))
    panel_b <- data.frame(phi_IDLE = plotted_phi, proteome_frac,
                          check.names = FALSE)[in_x, ]
    write_source(panel_b, fig_dir, "3B", "Figure_3B.csv")
    
    # Panel C: biomass composition (mass fractions, as plotted)
    biomass_frac <- biomass / rowSums(biomass)
    colnames(biomass_frac) <- sub("^c\\.", "", colnames(biomass_frac))
    panel_c <- data.frame(phi_IDLE = plotted_phi, biomass_frac,
                          check.names = FALSE)[in_x, ]
    write_source(panel_c, fig_dir, "3C", "Figure_3C.csv")
    
    writeLines(c(
      "Source data for Figure 3",
      "",
      paste0("Simulations: GBA model ", modelname, " (reversible kinetics, with Q sector)."),
      "Only converged simulations within the plotted range are included.",
      "",
      "3A/Figure_3A_model.csv",
      "  phi_IDLE: proteome fraction of the idle protein",
      "  relative_growth_rate: growth rate divided by the maximum growth rate (mu/mu*)",
      "",
      "3A/Figure_3A_Scott_2010.csv",
      "  Experimental data from Scott et al. (2010), Science 330:1099-1102 (gray triangles).",
      "  phi_IDLE: proteome fraction of the overexpressed protein (beta-galactosidase)",
      "  relative_growth_rate: growth rate divided by the growth rate without",
      "    overexpression in the same condition",
      "",
      "3B/Figure_3B.csv",
      "  Proteome composition: one column per protein; values are mass fractions",
      "  of the proteome (each row sums to 1).",
      "",
      "3C/Figure_3C.csv",
      "  Biomass composition: one column per biomass component; values are mass",
      "  fractions of the biomass (each row sums to 1)."
    ), here("figures", "source_data_figures", fig_dir, "README.txt"))
  }
}
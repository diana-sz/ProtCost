#### Plots for Dekel & Alon simulations ####
# Author: Diana Szeliova

library(RColorBrewer)
library(here)

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

predict.parameters <- 0
plotted_phis <- c(0.04, 0.3)
phi_xlim <- 0.5
cex_lab <- 0.8
cex_axis <- 1.1
ylim <- c(-3, 1)
cols_lines = c("grey15", "grey70")  # dashed lines

# dekel experimental data
dekel_conc <- c(0.0001, 0.001, 0.01, 0.1, 0.312, 0.6, 1, 3.875, 8.273, 12.09, 16.28)
dekel_rel_mu <- c(-4.79, -4.87, -6.11, -1.31, 0.93, 7.81, 9.66, 10.21, 13.76, 10.98, 12.37)
dekel_rel_mu <- 1+dekel_rel_mu/100


# ---- Plotting ----
plot_composition <- function(proteome, target_phi, colors = NULL,
                             main = NA,
                             ylab = NA,
                             xlab = NA,
                             legend = TRUE,
                             cex_lab = 1,
                             cex_axis = 1,
                             xaxis = TRUE,
                             yaxis = TRUE,
                             xlim = c(0,100),
                             ylim = c(0,1)) {
  
  # Cumulative sum by row (used for stacking)
  proteome <- proteome/rowSums(proteome)
  y_cum <- t(apply(proteome[ncol(proteome):1], 1, cumsum))
  
  # Base plot (empty)
  plot(NA, xlim = xlim, ylim =ylim, xlab = xlab, ylab = ylab, main = main,
       axes = FALSE, xaxs = "i", yaxs = "i",
       cex.lab = cex_lab, log = "x")
  if (xaxis) axis(1, cex.axis = cex_axis)
  if (yaxis) axis(2, las = 1, at = seq(0,1,0.2), cex.axis = cex_axis)
  
  box()
  
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
  
  # Add legend
  if(legend){
    leg_text <- gsub("c\\.|p\\.", "", colnames(proteome))
    par(xpd = NA)
    legend(xlim[1]*0.22, -0.25,
           legend = leg_text, fill = rev(colors), bty = "n", cex = cex_lab*1.1, ncol = 4)
    par(xpd = FALSE)
    
  }
}

modelname <- "M9_dekel_rev"
opt_data <- read.csv(here("data", paste0(modelname, ".csv")))
opt_data <- subset(opt_data, convergence == 4 & phi <= phi_xlim)

x_C2_vals <- sort(unique(opt_data$x_C2)) 
cols <- setNames(viridisLite::magma(length(x_C2_vals), end = 0.9), x_C2_vals)
prot_vals <- sort(unique(opt_data$convergence)) 
shapes <- setNames((20+length(prot_vals)-1):20, prot_vals)
xlim <- range(c(x_C2_vals, max(dekel_conc)))

opt_data$mu_ref <- opt_data$mu / max(opt_data[opt_data$phi == 0, "mu"])

# fit cost for the lowest x_C2 concentration
ref <- opt_data[opt_data$x_C2 == min(opt_data$x_C2),] 
ref_mu <- ref$mu 
fit <- lm(mu ~ phi, data=ref)

for (format in c("tif", "png")) {
  open_figure("Figure_4", 24, 8, format)
  
  # Panels 1 and 2 take full height, panel 3 is split into 4 subplots
  layout_matrix <- matrix(c(1,1,2,2,3,4), nrow = 2, byrow = FALSE)
  layout(layout_matrix, widths = c(0.8,0.8,0.45,0.45), heights = c(1,1))
  
  # --- Panel 1: Growth rate vs phi ---
  par(mar = c(4,3.5,2.5,1), oma = c(0.5,0.5,0,0.5))
  plot(NA,
       ylab = NA,
       xlab = NA,
       axes = FALSE,
       ylim = c(0, max(opt_data$mu_ref)),
       xlim = c(0, phi_xlim))
  # one line per L_ext concentration
  for (x2 in x_C2_vals) {
    one_x2 <- opt_data[opt_data$x_C2 == x2, ]
    one_x2 <- one_x2[order(one_x2$phi), ]
    lines(one_x2$phi, one_x2$mu_ref, col = cols[as.character(x2)], lwd = 3)
  }
  box()
  axis(1, cex.axis=cex_axis) 
  axis(2, las = 1, cex.axis=cex_axis) 
  mtext(expression("LAC proteome fraction " *  italic("\u03D5")["LAC"]), side=1, cex = cex_lab, line = 2.8)
  mtext(expression("Growth rate relative to " * italic("\u03D5")["LAC"] * "=0"), side=2, cex = cex_lab, line = 2.4)
  
  abline(v=plotted_phis, col = cols_lines, lty=2)
  legend("bottomleft", legend = x_C2_vals, col = cols, lwd = 3, title = expression("L"["ext"]))
  add_panel_label(LETTERS[1], cex = cex_lab*1.9)
  
  # --- Panel 2: Relative growth rate vs x_C2 for different phi ---
  plot(NA,
       ylab = NA, xlab = NA,
       log = "x",
       xlim = xlim,
       axes = FALSE,
       ylim = c(min(opt_data[opt_data$phi %in% plotted_phis, "mu_ref"]), 1.15))
  box()
  
  for(pphi in seq_along(plotted_phis)){
    one_phi <- opt_data[opt_data$phi == plotted_phis[pphi], ]
    point_cols <- cols[as.character(one_phi$x_C2)]
    lines(mu_ref ~ x_C2, data = one_phi, col = cols_lines[pphi], lty = 2)
    points(mu_ref ~ x_C2, data = one_phi, col = point_cols, pch = 20, cex = 1.6)
    text(min(x_C2_vals), (one_phi[which.min(one_phi$x_C2), "mu_ref"]+0.05), 
         bquote(italic("\u03D5")["LAC"] * " = " * .(plotted_phis[pphi])), cex = 0.9, adj=0)
  }
  
  points(dekel_conc, dekel_rel_mu, col = "grey50", pch = 24, cex = 1)
  
  axis(1, cex.axis=cex_axis) 
  axis(2, las = 1, cex.axis=cex_axis) 
  mtext(expression("Growth rate relative to " * italic("\u03D5")["LAC"] * "=0"), side=2, cex = cex_lab, line = 2.4)
  mtext(expression("L"["ext"] * " concentration"), side = 1, line = 2.8, cex = cex_lab)
  add_panel_label(LETTERS[2], cex = cex_lab*1.9)
  
  # --- Panel 3: 4 subplots for cost/benefit ---
  xlim_c <- c(xlim[1], max(x_C2_vals))
  
  cex_axis_small <- cex_axis*0.85
  for (p in seq_along(plotted_phis)) {    
    phi <- plotted_phis[p]
    biomass <- opt_data[opt_data$phi == phi, grep("c\\.", colnames(opt_data))]
    biomass <- biomass[, -grep("x_", colnames(biomass))]
    
    concs <- opt_data[opt_data$phi == phi, "x_C2"]
    colors2 <- rev(brewer.pal(ncol(biomass), "RdBu"))
    plot_legend <- ifelse(p == 2, TRUE, FALSE)
    draw_axis <- ifelse(p == 2, TRUE, FALSE)
    
    upper <- ifelse(p == 1, 2.5, 0)
    par(mar = c(4-upper, 3.5, upper, 0.5))
    plot_composition(biomass, concs, colors2, main="", 
                     ylab = NA,
                     xlab = NA,
                     cex_lab=cex_lab,
                     cex_axis=cex_axis_small,
                     legend=plot_legend,
                     xaxis = draw_axis,
                     xlim = xlim_c)
    mtext(bquote(italic("\u03D5")["LAC"] * " = " * .(phi)), side=3, cex = 0.6, line = 0.12, adj = 1)
    
    if(p == 1){
      add_panel_label(LETTERS[3], cex = cex_lab*1.9)
    }
    
  }
  mtext("Biomass composition", side = 2, line = 2.4, cex = cex_lab, outer = FALSE, adj = -1.7) 
  
  dev.off()
}


# ---- Source data ----
source_dir <- here("figures", "source_data_figures", "Figure 4")

write_source <- function(df, panel, file_name) {
  out_dir <- file.path(source_dir, panel)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(out_dir, file_name), row.names = FALSE)
}

# Panel A: relative growth rate vs LAC proteome fraction, for each L_ext
panel_a <- data.frame(L_ext = opt_data$x_C2,
                      phi_LAC = opt_data$phi,
                      relative_growth_rate = opt_data$mu_ref)
panel_a <- panel_a[order(panel_a$L_ext, panel_a$phi_LAC), ]
write_source(panel_a, "4A", "Figure_4A.csv")

# Panel B: simulations at the plotted phi values, and experimental data
panel_b_sim <- panel_a[panel_a$phi_LAC %in% plotted_phis,
                       c("phi_LAC", "L_ext", "relative_growth_rate")]
panel_b_sim <- panel_b_sim[order(panel_b_sim$phi_LAC, panel_b_sim$L_ext), ]
write_source(panel_b_sim, "4B", "Figure_4B_simulations.csv")

panel_b_exp <- data.frame(L_ext = dekel_conc,
                          relative_growth_rate = dekel_rel_mu)
write_source(panel_b_exp, "4B", "Figure_4B_Dekel_Alon_2005.csv")

# Panel C: biomass composition (mass fractions, as plotted)
panel_c <- do.call(rbind, lapply(plotted_phis, function(phi) {
  rows <- opt_data$phi == phi
  biomass <- opt_data[rows, grep("c\\.", colnames(opt_data))]
  biomass <- biomass[, -grep("x_", colnames(biomass))]
  biomass <- biomass / rowSums(biomass)
  colnames(biomass) <- sub("^c\\.", "", colnames(biomass))
  out <- data.frame(phi_LAC = phi, L_ext = opt_data$x_C2[rows], biomass,
                    check.names = FALSE)
  out[order(out$L_ext), ]
}))
write_source(panel_c, "4C", "Figure_4C.csv")

# README describing the files and columns
writeLines(c(
  "Source data for Figure 4",
  "",
  "Simulations: GBA model M9_dekel_rev (reversible kinetics).",
  "Only converged simulations with phi_LAC <= 0.5 are included, as in the figure.",
  "",
  "4A/Figure_4A.csv",
  "  L_ext: external lactose concentration",
  "  phi_LAC: proteome fraction of the LAC protein",
  "  relative_growth_rate: growth rate divided by the maximum growth rate at phi_LAC = 0",
  "",
  "4B/Figure_4B_simulations.csv",
  "  Simulated relative growth rate at phi_LAC = 0.04 and 0.3 (dashed lines).",
  "  Columns as in 4A.",
  "",
  "4B/Figure_4B_Dekel_Alon_2005.csv",
  "  Experimental data from Dekel and Alon (2005), Nature 436:588-592 (gray triangles).",
  "  relative_growth_rate: 1 + (reported relative growth change in %) / 100",
  "",
  "4C/Figure_4C.csv",
  "  Simulated biomass composition at phi_LAC = 0.04 and 0.3.",
  "  One column per biomass component; values are mass fractions (each row sums to 1)."
), file.path(source_dir, "README.txt"))
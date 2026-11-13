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

models <- c("M9_dekel_rev", "M10_dekel_efflux_rev")
predict.parameters <- 0
# Names of the figure files in the manuscript; other outputs keep their model name
figure_names <- c(M9_dekel_rev = "Appendix_Figure_S11_no_pump",
                  M10_dekel_efflux_rev = "Appendix_Figure_S11_pump")
plotted_phis <- c(0.04, 0.3)
phi_xlim <- 0.5
cex_all <- 0.9
ylim <- c(-3, 1)
cex_lab <- 0.85
cex_axis <- 1.1
source(here("code", "uni_colors.R"))

# dekel experimental data
dekel_conc <- c(0.0001, 0.001, 0.01, 0.1, 0.312, 0.6, 1, 3.875, 8.273, 12.09, 16.28)
dekel_rel_mu <- c(-4.79, -4.87, -6.11, -1.31, 0.93, 7.81, 9.66, 10.21, 13.76, 10.98, 12.37)
dekel_rel_mu <- 1+dekel_rel_mu/100


# ---- Plotting ----
plot_composition <- function(proteome, target_phi,
                             colors,
                             xlim,
                             ylim = c(0,1),
                             main = "",
                             ylab = NA,
                             xlab = NA,
                             legend = TRUE,
                             cex_lab = 1,
                             cex_axis = 1,
                             xaxis = TRUE,
                             yaxis = TRUE) {
  
  # Cumulative sum by row (used for stacking)
  proteome <- proteome/rowSums(proteome)
  y_cum <- t(apply(proteome[ncol(proteome):1], 1, cumsum))
  
  # Base plot (empty)
  plot(NA, xlim = xlim, ylim =ylim, xlab = xlab, ylab = ylab, main = main,
       xaxs = "i", yaxs = "i",
       axes = FALSE, cex.lab = cex_lab, log = "x")
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
    legend(xlim[1], -0.25,
           legend = leg_text, fill = rev(colors), bty = "n", cex = cex_lab*1.15, ncol = 4)
    par(xpd = FALSE)
    
  }
}


for(model in models){
  is.reversible <- ifelse(grepl("rev", model), 1, 0)
  modelname <- gsub("_rev", "", model)
  
  source(here("code", "initialize_model.R"))
  
  opt_data <- read.csv(here("data", paste0(modelname, ".csv")))
  opt_data <- subset(opt_data, convergence == 4 & phi <= phi_xlim)
  
  row <- 1
  rho_cond <- rho_cond[1]
  a_cond  <- a_cond[,1, drop=FALSE]
  
  # ## calculate cost contributions
  res <- vector("list", 0)
  M_colSums <- colSums(M)
  
  for (phi in plotted_phis){ #one_prot$phi
    one_phi <- opt_data[opt_data$phi == phi,]
    
    for (x_C2 in unique(one_phi$x_C2)) {
      one_xC2 <- one_phi[one_phi$x_C2 == x_C2,]
      
      a_cond[1,] <- one_xC2$x_C[1]
      a_cond[3,] <- x_C2
      
      vs <- one_xC2[row, grep("v\\.", colnames(one_xC2))]
      fs <- one_xC2[row, grep("f\\.", colnames(one_xC2))]
      
      cint <- one_xC2[row, grep("c\\.", colnames(one_xC2))]
      cint <- cint[, (1+length(a_cond)):length(cint)]
      fint <- cint/rho_cond
      taus <- tau(a_cond, t(cint))
      dtaus <- dtau(a_cond, t(cint))
      growth_rate <- one_xC2[row, "mu"]
      
      for(j in 1:length(vs)){
        Mjp <- M[nrow(M),j]  # protein row
        catalytic_cost <- growth_rate*taus[j]
        saturation_value <- unlist(vs) %*% (dtaus %*% M[, j])
        crowding_value <- M_colSums[j] * unlist(vs) %*% dtaus %*% unlist(fint)
        
        res[[length(res) + 1]] <- data.frame(
          growth_rate = one_xC2$mu,
          x_C2 = x_C2,
          phi = phi,
          variable = colnames(vs)[j],
          protein_benefit = as.numeric(Mjp),
          catalytic_cost = -as.numeric(catalytic_cost),
          saturation_value = -as.numeric(saturation_value),
          crowding_value = as.numeric(crowding_value),
          sum = as.numeric(Mjp - catalytic_cost - saturation_value + crowding_value)
        )
      }
    }
    results <- do.call(rbind, res)
  }
  
  plots <- c("sum", "catalytic_cost", "saturation_value", "crowding_value")
  plot_names <- c("marginal value", "catalytic cost",
                  "saturation value", "crowding value")
  colors <- c("black", uni_red, uni_green, uni_lila)
  x_C2_vals <- sort(unique(opt_data$x_C2)) 
  cols <- setNames(viridisLite::magma(length(x_C2_vals), end = 0.9), x_C2_vals)
  xlim <- range(c(x_C2_vals, max(dekel_conc)))
  
  opt_data$mu_ref <- opt_data$mu / max(opt_data[opt_data$phi == 0, "mu"])
  
  file_name <- model
  if (file_name %in% names(figure_names)) file_name <- figure_names[[file_name]]
  for (format in c("tif", "png")) {
    open_figure(file_name, 26, 8, format)
    
    # Panels 1 and 2 take full height, panel 3 is split into 4 subplots
    layout_matrix <- matrix(c(1,1,2,2,3,4,5,6), nrow = 2, byrow = FALSE)
    layout(layout_matrix, widths = c(0.8,0.8,0.45,0.45), heights = c(1,1))
    
    # --- Panel 1: Growth rate vs phi ---
    par(mar = c(4,4,3,1), oma = c(1,0.5,0,0))
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
    mtext(expression("LAC proteome fraction " * italic("\u03D5")["LAC"]), side=1, cex = cex_lab, line = 2.8)
    mtext(expression("Growth rate relative to " * italic("\u03D5")["LAC"] * "=0"), side=2, cex = cex_lab, line = 2.4)
    cols_lines = c("grey15", "grey60")
    abline(v=plotted_phis, col = cols_lines, lty=2)
    legend("bottomleft", legend = x_C2_vals, col = cols, lwd = 3, title = expression("L"["ext"]))
    add_panel_label(LETTERS[1], cex = cex_lab*1.9)
    
    
    # --- Panel 2: Relative growth rate vs x_C2 for different phi ---
    plot(NA,
         ylab = NA, xlab = NA,
         log = "x",
         xlim = xlim,
         axes = FALSE,
         ylim = c(0.65, 1.15)) #max(opt_data$mu_ref)))
    box()
    
    for(pphi in seq_along(plotted_phis)){
      one_phi <- opt_data[opt_data$phi == plotted_phis[pphi], ]
      point_cols <- cols[as.character(one_phi$x_C2)]
      lines(mu_ref ~ x_C2, data = one_phi, col = cols_lines[pphi], lty = 2)
      points(mu_ref ~ x_C2, data = one_phi, col = point_cols, pch = 20, cex = 1.6)
      text(min(x_C2_vals), (one_phi[which.min(one_phi$x_C2), "mu_ref"]+0.055), 
           bquote(italic("\u03D5")["LAC"] * " = " * .(plotted_phis[pphi])), cex = cex_all*1.1, adj=0)
    }
    
    points(dekel_conc, dekel_rel_mu, col = "grey50", pch = 24, cex = 1)
    
    axis(1, cex.axis=cex_axis) 
    axis(2, las = 1, cex.axis=cex_axis) 
    mtext(expression("Growth rate relative to " * italic("\u03D5")["LAC"] * "=0"), side=2, cex = cex_lab, line = 2.4)
    mtext(expression("L"["ext"] * " concentration [g/L]"), side = 1, line = 2.8, cex = cex_lab)
    add_panel_label(LETTERS[2], cex = cex_lab*1.9)
    
    # --- Panel 3: 4 subplots for cost/benefit ---
    xlim_c <- c(xlim[1], max(x_C2_vals))
    
    cex_axis_small <- cex_axis * 0.85
    for (p in seq_along(plotted_phis)) {
      phi <- plotted_phis[p]
      one_phi <- results[results$phi == phi & results$variable == "v.LAC", ]
      
      ls <- ifelse(p %in% c(2), 1.4, 3.7)
      par(mar = c(2.5, ls, 3, 4.2-ls))
      
      plot(NA, xlim = xlim_c, ylim = ylim, log = "x", axes = FALSE, xlab = NA, ylab = NA)
      box()
      abline(h = 0, col = "grey80")
      
      for (i in length(plots):1) {
        lines(one_phi$x_C2, one_phi[[plots[i]]], col = colors[i], lty = i, lwd = 1.5)
      }
      mtext(bquote(italic("\u03D5")["LAC"] * " = " * .(phi)), side=3, cex = 0.75, line = 0.2)
      if (p == 1) add_panel_label(LETTERS[3], cex = cex_lab*1.9)
      
      if(p == 1){
        axis(2, las = 1, at = seq(ylim[1], ylim[2], 1), cex.axis=cex_axis_small) 
        mtext("Cost / Benefit", side = 2, line = 2.4, cex = cex_lab)
        par(xpd = NA) 
        legend(xlim[1]*0.9, ylim[1]*0.95,
               legend = rev(plot_names), col = rev(colors), cex = cex_lab*1.15, lwd = 1.5,
               lty = length(plots):1, bty = "n", horiz=FALSE, ncol=2)
        par(xpd = FALSE)
      }
      
      biomass <- opt_data[opt_data$phi == phi, grep("c\\.", colnames(opt_data))]
      biomass <- biomass[, -grep("x_", colnames(biomass))]
      
      concs <- opt_data[opt_data$phi == phi, "x_C2"]
      colors2 <- rev(brewer.pal(ncol(biomass), "RdBu"))
      plot_legend <- ifelse(p == 1, TRUE, FALSE)
      draw_axis <- ifelse(p == 1, TRUE, FALSE)
      
      par(mar = c(4, ls, 0, 4-ls))
      plot_composition(biomass, concs, colors2, main="", 
                       ylab = NA,
                       xlab = NA,
                       cex_lab=cex_lab,
                       legend=plot_legend,
                       yaxis = draw_axis,
                       xlim = xlim_c)
      if(p == 1){
        mtext("Biomass composition", side = 2, line = 2.4, cex = cex_lab) 
      }
    }
    dev.off()
  }
}
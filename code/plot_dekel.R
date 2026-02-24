#### Plots for Dekel & Alon simulations ####
# Author: Diana Szeliova

library(RColorBrewer)
library(here)

models <- c("M9_dekel", "M10_dekel_efflux")
predict.parameters <- 0
plotted_phis <- c(0.04, 0.5)
phi_xlim <- 0.5
cex_all <- 0.9
ylim <- c(-3, 1)
cex_lab <- 0.8
cex_axis <- 1.1
letters_line <- 1.25
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
  y_cum <- t(apply(proteome[m:1], 1, cumsum))
  
  # Base plot (empty)
  plot(NA, xlim = xlim, ylim =ylim, xlab = xlab, ylab = ylab, main = main,
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
    legend(xlim[1], -0.33,  #"bottomleft", 
           legend = leg_text, fill = rev(colors), bty = "n", cex = cex_lab, ncol = 4)
    par(xpd = FALSE)
    
  }
}


for(is.reversible in c(1,0)){
  for(modelname in models){
    
    source(here("code", "initialize_model.R"))
    
    opt_data <- read.csv(here("data", paste0(modelname, ".csv")))
    opt_data <- subset(data, convergence == 4 & phi <= phi_xlim)
    
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
    cols <- setNames(rev(brewer.pal(length(x_C2_vals), "Paired")), x_C2_vals) 
    prot_vals <- sort(unique(opt_data$convergence)) 
    shapes <- setNames((20+length(prot_vals)-1):20, prot_vals)
    xlim <- range(c(x_C2_vals, max(dekel_conc)))
    
    opt_data$mu_ref <- opt_data$mu / max(opt_data[opt_data$phi == 0, "mu"])
    
    png(here("figures", paste0(modelname, is.reversible, ".png")),
        type="cairo", units="cm",
        width=26, height=8, res=300)
    
    # Panels 1 and 2 take full height, panel 3 is split into 4 subplots
    layout_matrix <- matrix(c(1,1,2,2,3,4,5,6), nrow = 2, byrow = FALSE)
    layout(layout_matrix, widths = c(0.8,0.8,0.45,0.45), heights = c(1,1))
    
    # --- Panel 1: Growth rate vs phi ---
    par(mar = c(4,4,3,1), oma = c(1,0.5,0,0))
    point_cols <- cols[as.character(opt_data$x_C2)] 
    plot(mu_ref ~ phi, data = opt_data, 
         ylab = NA,
         xlab = NA,
         axes = FALSE,
         ylim = c(0, max(mu_ref)),
         xlim = c(0, phi_xlim),
         cex = 1,
         col = point_cols, pch = 20)
    box()
    axis(1, cex.axis=cex_axis) 
    axis(2, las = 1, cex.axis=cex_axis) 
    mtext(expression("LAC proteome fraction (" * Phi * ")"), side=1, cex = cex_lab, line = 2.8)
    mtext(expression("Growth rate relative to " * Phi["LAC"] * "=0"), side=2, cex = cex_lab, line = 2.4)
    cols_lines = c("grey15", "grey60")
    abline(v=plotted_phis, col = cols_lines, lty=2)
    legend("bottomleft", legend = x_C2_vals, col = cols, pch = 16, title = expression("L"["ext"]))
    mtext(
      paste0("(", letters[1], ")"),
      side = 3,
      adj = 0.5,
      line = letters_line,
      cex = cex_lab,
      font = 2
    )
    
    
    # --- Panel 2: Relative growth rate vs x_C2 for different phi ---
    plot(NA,
         ylab = NA, xlab = NA,
         log = "x",
         xlim = xlim,
         axes = FALSE,
         ylim = c(min(opt_data[opt_data$phi %in% plotted_phis, "mu_ref"]), 1.15)) #max(opt_data$mu_ref)))
    box()
    
    for(pphi in seq_along(plotted_phis)){
      one_phi <- opt_data[opt_data$phi == plotted_phis[pphi], ]
      point_cols <- cols[as.character(one_phi$x_C2)]
      lines(mu_ref ~ x_C2, data = one_phi, col = cols_lines[pphi], lty = 2)
      points(mu_ref ~ x_C2, data = one_phi, col = point_cols, pch = 20, cex = 1.3)
      text(min(x_C2_vals), (one_phi[which.min(one_phi$x_C2), "mu_ref"]+0.055), 
           bquote(Phi * " = " * .(plotted_phis[pphi])), cex = cex_all, adj=0)
    }
    
    points(dekel_conc, dekel_rel_mu, col = "grey50", pch = 24, cex = 1)
    
    axis(1, cex.axis=cex_axis) 
    axis(2, las = 1, cex.axis=cex_axis) 
    mtext(expression("Growth rate relative to " * Phi["LAC"] * "=0"), side=2, cex = cex_lab, line = 2.4)
    mtext(expression("L"["ext"] * " concentration"), side = 1, line = 2.8, cex = cex_lab)
    mtext(
      paste0("(", letters[2], ")"),
      side = 3,
      adj = 0.5,
      line = letters_line,
      cex = cex_lab,
      font = 2
    )
    
    # --- Panel 3: 4 subplots for cost/benefit ---
    xlim[2] <- max(x_C2_vals)
    mtext(
      paste0("(", letters[3], ")"),
      side = 3,
      adj = 1.89,
      line = letters_line,
      cex = cex_lab,
      font = 2
    )
    
    cex_axis_small <- cex_axis * 0.85
    for (p in seq_along(plotted_phis)) {
      phi <- plotted_phis[p]
      one_phi <- results[results$phi == phi & results$variable == "v.LAC", ]
      
      ls <- ifelse(p %in% c(2), 1.4, 3.7)
      par(mar = c(2.5, ls, 3, 4.2-ls))
      
      plot(NA, xlim = xlim, ylim = ylim, log = "x", axes = FALSE, xlab = NA, ylab = NA)
      box()
      abline(h = 0, col = "grey80")
      
      for (i in length(plots):1) {
        lines(one_phi$x_C2, one_phi[[plots[i]]], col = colors[i], lty = i, lwd = 1.5)
      }
      mtext(bquote(Phi * " = " * .(phi)), side=3, cex = 0.7, line = 0.3)
      
      if(p == 1){
        axis(2, las = 1, at = seq(ylim[1], ylim[2], 1), cex.axis=cex_axis_small) 
        mtext("Cost / Benefit", side = 2, line = 2.4, cex = cex_lab)
        par(xpd = NA) 
        legend(xlim[1], ylim[1], #inset = c(-0.05,-0.05), 
               legend = rev(plot_names), col = rev(colors), cex = cex_lab, lwd = 1.5,
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
                       xlim = xlim)
      if(p == 1){
        mtext("Biomass composition", side = 2, line = 2.4, cex = cex_lab) 
      }
    }
    dev.off()
  }
}
 


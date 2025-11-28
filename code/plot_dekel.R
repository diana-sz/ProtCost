library(RColorBrewer)
library(here)
library(readODS)

directory <- paste0(here(), "/code")
setwd(directory)

predict.parameters <- 0
plotted_phis <- c(0.04, 0.4)
source("uni_colors.R")


# ---- Plotting ----
plot_composition <- function(proteome, target_phi, colors = NULL,
                             main = "Proteome Composition",
                             ylab = "Composition",
                             xlab = "Proteome fr
                             action",
                             legend = TRUE,
                             cex_lab = 1,
                             xaxis = TRUE,
                             yaxis = TRUE,
                             xlim = c(0,100), ylim = c(0,1)) {
  # Dimensions
  n <- nrow(proteome)
  m <- ncol(proteome)
  
  # Cumulative sum by row (used for stacking)
  proteome <- proteome/rowSums(proteome)
  y_cum <- t(apply(proteome[ncol(proteome):1], 1, cumsum))
  
  # Base plot (empty)
  plot(NA, xlim = xlim, ylim =ylim, xlab = xlab, ylab = ylab, main = main,
       axes = FALSE,
       cex.lab = cex_lab, log = "x")
  if(yaxis){axis(2)}
  axis(1)
  box()
  
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
    par(xpd = NA)
    legend(xlim[1], -0.33,  #"bottomleft", 
           legend = leg_text, fill = rev(colors), bty = "n", cex = 1.05, ncol = 4)
    par(xpd = FALSE)
    
  }
}



for(is.reversible in c(1,0)){
  modelname <- "A8dekel"
  source("initialize_model.R")
  
  opt_data <- read.csv(paste0("../data/", modelname, ".csv"), row.names = 1)
  opt_data <- opt_data[opt_data$convergence == 4, ]
  
  row <- 1
  rho_cond <- rho_cond[1]
  a_cond  <- a_cond[,1, drop=FALSE]
  
  # ## calculate cost contributions
  results <- data.frame(
    growth_rate = numeric(),
    x_C2 = numeric(),
    phi = numeric(),
    variable = character(),
    protein_benefit = numeric(),
    local_cost = numeric(),
    local_benefit = numeric(),
    transport_benefit = numeric(),
    sum = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (phi in plotted_phis){ #one_prot$phi
    one_phi <- opt_data[opt_data$phi == phi,]
    
    for (x_C2 in unique(one_phi$x_C2)) {
      one_xC2 <- one_phi[one_phi$x_C2 == x_C2,]
      
      a_cond[1,] <- one_xC2$x_C[1]
      a_cond[2,] <- x_C2
      
      vs <- one_xC2[row, grep("v\\.", colnames(one_xC2))]
      fs <- one_xC2[row, grep("f\\.", colnames(one_xC2))]
      
      cint <- one_xC2[row, grep("c\\.", colnames(one_xC2))]
      cint <- cint[, (1+length(a_cond)):length(cint)]
      fint <- cint/rho_cond
      taus <- tau(a_cond, t(cint))
      dtaus <- dtau(a_cond, t(cint))
      growth_rate <- one_xC2[row, "mu"]
      
      for(j in 1:length(vs)){
        Mjp <- M["p",j]
        local_cost <- growth_rate*taus[j]
        local_benefit <- unlist(vs) %*% (dtaus %*% M[, j])
        transport_benefit <- colSums(M)[j] * unlist(vs) %*% dtaus %*% unlist(fint)
        
        results <- rbind(
          results,
          data.frame(
            growth_rate = one_xC2$mu,
            x_C2 = x_C2,
            phi = phi,
            variable = colnames(vs)[j],
            protein_benefit = as.numeric(Mjp),
            local_cost = -as.numeric(local_cost),
            local_benefit = -as.numeric(local_benefit),
            transport_benefit = as.numeric(transport_benefit),
            sum = as.numeric(Mjp - local_cost - local_benefit + transport_benefit)
          )
        )
      }
    }
  }
  

  plots <- c("sum", "local_cost", "local_benefit", "transport_benefit")
  plot_names <- c("sum", "(marginal) protein investment",
                  "(marginal) local protein value", "biomass production cost")
  colors <- c(uni_red, uni_green, uni_lila, uni_blue)
  x_C2_vals <- sort(unique(opt_data$x_C2)) 
  cols <- setNames(rev(brewer.pal(length(x_C2_vals), "Paired")), x_C2_vals) 
  prot_vals <- sort(unique(opt_data$convergence)) 
  shapes <- setNames((20+length(prot_vals)-1):20, prot_vals)
  
  # dekel data
  dekel_conc <- c(0.0001, 0.001, 0.01, 0.1, 0.312, 0.6, 1, 3.875, 8.273, 12.09, 16.28)
  dekel_rel_mu <- c(-4.79, -4.87, -6.11, -1.31, 0.93, 7.81, 9.66, 10.21, 13.76, 10.98, 12.37)
  dekel_rel_mu <- 1+dekel_rel_mu/100
  
  
  xlim <- range(c(x_C2_vals, max(dekel_conc)))
  ylim <- c(-3, 1)
  cex_lab <- 1.05
  cex_axis <- 1.15
  
  opt_data$mu_ref <- opt_data$mu / max(opt_data[opt_data$phi == 0, "mu"])
  
  # fit cost for the lowest x_C2 concentration
  ref <- opt_data[opt_data$x_C2 == min(opt_data$x_C2),] 
  ref_mu <- ref$mu 
  fit <- lm(mu ~ phi, data=ref)
  
  png(paste0("../figures/", modelname, is.reversible, ".png"), 
      type="cairo", units="cm",
      width=26, height=8, res=300)
  
  # Panels 1 and 2 take full height, panel 3 is split into 4 subplots
  layout_matrix <- matrix(c(1,1,2,2,3,4,5,6), nrow = 2, byrow = FALSE)
  layout(layout_matrix, widths = c(0.8,0.8,0.45,0.45), heights = c(1,1))
  
  # --- Panel 1: Growth rate vs phi ---
  par(mar = c(4,4,2,1), oma = c(1,0.5,0,0))
  point_cols <- cols[as.character(opt_data$x_C2)] 
  #point_shapes <- shapes[as.character(opt_data$convergence)]
  plot(mu ~ phi, data = opt_data, 
       ylab = NA,
       xlab = NA,
       axes = FALSE,
       ylim = c(0, max(mu)),
       cex = 1,
       col = point_cols, pch = 20)
  box()
  axis(1, cex.axis=cex_axis) 
  axis(2, las = 1, cex.axis=cex_axis) 
  mtext(expression("LAC proteome fraction (" * Phi * ")"), side=1, cex = cex_lab, line = 2.8)
  mtext(expression("Growth rate [h"^-1 * "]"), side=2, cex = cex_lab, line = 2.4)
  legend("topright", legend = x_C2_vals, col = cols, pch = 16, title = "xL")
  
  abline(fit)
  cols_lines = c("grey20", "grey55")
  abline(v=plotted_phis, col = cols_lines, lty=2)
  
  # --- Panel 2: Relative growth rate vs x_C2 for different phi ---
  #par(mar = c(4,4,2,1))
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
    #point_shapes <- shapes[as.character(one_phi$convergence)]
    lines(mu_ref ~ x_C2, data = one_phi, col = cols_lines[pphi], lty = 2)
    points(mu_ref ~ x_C2, data = one_phi, col = point_cols, pch = 20, cex = 1.3)
    text(min(x_C2_vals), (one_phi[which.min(one_phi$x_C2), "mu_ref"]+0.055), 
         bquote(Phi * " = " * .(plotted_phis[pphi])), cex = 1.05, adj=0)
  }
  
  
  points(dekel_conc, dekel_rel_mu, col = "grey50", pch = 24, cex = 1.1)
  
  axis(1, cex.axis=cex_axis) 
  axis(2, las = 1, cex.axis=cex_axis) 
  #mtext("xC2 concentration", side=1, cex = cex_lab, line = 2.5)
  mtext(expression("Growth rate relative to " * Phi["tC2"] * "=0"), side=2, cex = cex_lab, line = 2.4)
  mtext(expression("L"["ext"] * " concentration"), side = 1, line = 2.8, cex = cex_lab)
  
  
  # --- Panel 3: 4 subplots for cost/benefit ---
  for (p in seq_along(plotted_phis)) {
    phi <- plotted_phis[p]
    one_phi <- results[results$phi == phi & results$variable == "v.LAC", ]
    
    ls <- ifelse(p %in% c(2), 1, 3)
    par(mar = c(2.5, ls, 2, 4-ls))
    
    plot(NA, xlim = xlim, ylim = ylim, log = "x", axes = FALSE, xlab = NA, ylab = NA)
    box()
    abline(h = 0, col = "grey80")
    
    for (i in length(plots):1) {
      lines(one_phi$x_C2, one_phi[[plots[i]]], col = colors[i], lty = i, lwd = 1.5)
    }
    mtext(bquote(Phi * " = " * .(phi)), side=3, cex = cex_lab*0.85, line = 0.3)
    
    if (p %in% c(1,3)) { 
      axis(2, las = 1, at = seq(ylim[1], ylim[2], 1), cex.axis=cex_axis) 
      mtext("Cost / Benefit", side = 2, line = 2.2, cex = cex_lab*0.8)
      par(xpd = NA) 
      legend(xlim[1]-0.0008, ylim[1], #inset = c(-0.05,-0.05), 
             legend = rev(plot_names), col = rev(colors), cex = 0.85, lwd = 2,
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
                     xlim = c(min(concs), max(concs)))
    if(p == 1){
      mtext("Biomass composition", side = 2, line = 2.2, cex = cex_lab*0.8) 
      #mtext("LAC concentration", side = 1, line = 2.8, cex = cex_lab, adj = -70)
    }
  }
  
  dev.off()
  
}




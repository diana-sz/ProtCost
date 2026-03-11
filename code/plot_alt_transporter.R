library(here)
library(RColorBrewer)

cex_all <- 0.97
cex_points <- 0.95

source(here("code", "uni_colors.R"))

experimental <- data.frame(
  titrated_level = c(
    0,
    0.111111111,
    0.222222222,
    0.444444444,
    0.666666667,
    0.888888889,
    1.333333333,
    2.666666667
  ),
  mu_rel = c(
    0.46875,
    0.625,
    0.8125,
    0.9375,
    0.9375,
    1.0,
    1.03125,
    0.90625
  )
)

modelnames <- c("M9_alt_trans")
xlim <- c(0, 3)
ylim <- c(0, 1.1)

for(modelname in modelnames){
  results <- read.csv(here("data", paste0(modelname, ".csv")))
  results <- results[results$convergence == 4,]  # plot only simulations that converged
  
  png(here("figures", paste0(modelname,".png")),
      type="cairo", units="cm", res = 300,
      width=9, height=7.5)
  
  par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(3.7, 3.7, 0.5, 1))

  # Set factor levels
  opt_phi <- which.max(results$mu[results$kcat == min(results$kcat) & results$km == min(results$km)])
  results$rel_phi <- results$phi / results$phi[opt_phi]
  results$mu_rel <- results$mu / results$mu[opt_phi]
  results <- results[results$rel_phi <= 3,]
  
  # Set colors based on kcat
  kcat_vals <- sort(unique(results$kcat))
  palette_colors <- brewer.pal(length(kcat_vals), "Paired")
  color_map <- setNames(palette_colors, kcat_vals)
  point_colors <- color_map[as.character(results$kcat)]
  
  # Set shapes based on km 
  km_vals <- sort(unique(results$km))
  shape_list <- 19:(20 - length(km_vals))
  shape_map <- setNames(shape_list, km_vals)
  point_shapes <- shape_map[as.character(results$km)]
  
  # Create the plot
  plot(
    mu_rel ~ rel_phi, data = results,
    xlab = NA,
    ylab = NA,
    xlim = xlim, 
    ylim = ylim,
    pch = point_shapes,
    col = point_colors,
    cex = cex_points,
    axes = FALSE)

  mtext("Relative growth rate", side = 2, outer = FALSE, line = 2.6, cex = cex_all)
  mtext(expression("TS proteome fraction relative to optimum " * italic("\u03A6") * "/" * italic("\u03A6")^"\u204E"), 
                   side = 1, outer = FALSE, line = 2.5, cex = cex_all, adj = 0.75)
  
  box()
  axis(1, cex.axis = cex_all)
  axis(2, las=2, cex.axis = cex_all)
  
  # Add legend
  legend(
    0.9, 0, yjust = 0,
    legend = kcat_vals,
    col = palette_colors,
    pch = 15,
    cex = cex_points,
    bty = "n",
    title = expression(k[cat]^"TS2")
  )
  
  # Add legend
  legend(
    1.7, 0, yjust = 0,
    legend = km_vals,
    col = rep("black", length(unique(results$km))),
    pch = shape_list,
    cex = cex_points,
    bty = "n",
    title = expression(K[TS2]^"Sext")
  )

  dev.off()
  
}

png(here("figures", "experimental_transporter.png"),
    type="cairo", units="cm", res = 300,
    width=9, height=7.5)

# experimental data
#  PTS transporter system (glucose-specific subunit IIA) in Salmonella typhimurium

par(mfrow=c(1,1), oma = c(0,0,0,0))
par(mar = c(3.7, 3.7, 0.5, 1))

plot(
  mu_rel ~ titrated_level, data = experimental, 
  xlab = NA,
  ylab = NA,
  xlim = xlim,
  ylim = ylim,
  pch = 19,
  col = "black",
  axes = FALSE,
  cex = cex_points
)

x <- seq(0, 3, length.out = 50)
y <- 0.15/0.32 + 1.25 * (x / (0.5 + x))^1.2 - 0.22 * x
lines(x, y, col = "black")

box()
axis(1, cex.axis = cex_all)
axis(2, las=2, cex.axis = cex_all)
mtext("Relative growth rate", side = 2, outer = FALSE, line = 2.6, cex = cex_all)
mtext("Titrated level / WT level of transporter", side = 1, outer = FALSE, line = 2.5, cex = cex_all)

legend("bottomright", legend = "Maltose", col = "black", pch = 19, bty = "n", cex=cex_all)


dev.off()


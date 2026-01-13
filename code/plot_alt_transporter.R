library(here)
library(RColorBrewer)

setwd(here())

cex_all <- 1.05

source("code/uni_colors.R")

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



modelnames <- c("M9_alt_trans")#, "M9_alt_trans_rev")
xlim <- c(0, 3)
ylim <- c(0, 1.1)

for(modelname in modelnames){
  results <- read.csv(paste0("data/", modelname, ".csv"))
  results <- results[results$convergence == 4,]  # plot only simulations that converged
  
  png(paste0("figures/", modelname,".png"),
      type="cairo", units="cm", res = 300,
      width=17, height=7.5)
  
  par(mfrow=c(1,2), oma = c(0,1,0.5,0))
  par(mar = c(4, 2.5, 1, 1))

  # Set factor levels
  results$kcat <- results$kcat
  results$kcat <- factor(results$kcat, levels = unique(results$kcat))
  results$km <- factor(results$km, levels = unique(results$km))
  opt_phi <- which.max(results[results$kcat == 50.6 & results$km == 1, "mu"])
  results$rel_phi <- results$phi / results$phi[opt_phi]
  results$mu_rel <- results$mu / results$mu[opt_phi]
  results <- results[results$rel_phi <= 3,]
  
  # Set colors based on kcat
  palette_colors <- brewer.pal(length(unique(results$kcat)), "Paired")
  color_map <- setNames(palette_colors, unique(results$kcat))
  point_colors <- color_map[as.character(results$kcat)]
  
  # Set shapes based on km 
  shape_list <- 19:(18 - length(unique(results$km)) + 1)
  shape_map <- setNames(shape_list, unique(results$km))
  point_shapes <- shape_map[as.character(results$km)]
  
  # Create the plot
  plot(
    mu_rel ~ rel_phi, data = results,
    xlab = NA,
    ylab = NA,
    xlim = xlim, ylim = ylim,
    pch = point_shapes,
    col = point_colors,
    cex = cex_all, cex.lab = cex_all,
    axes = FALSE
    #main = "Simulations"
  )
  
  mtext("Relative growth rate", side = 2, outer = FALSE, line = 2.5, cex = cex_all)
  mtext("Fixed level / optimal level of TC", side = 1, outer = FALSE, line = 2.5, cex = cex_all)
  
  box()
  axis(1)
  axis(2, las=2)
  
  # Add legend
  legend(
    1, 0, yjust = 0,
    legend = unique(results$kcat),
    col = palette_colors,
    pch = 15,
    cex = 0.9,
    bty = "n",
    title = expression(k[cat]^"TC2")
  )
  
  # Add legend
  legend(
    2, 0, yjust = 0,
    legend = c(unique(results$km)),
    col = rep("black", length(unique(results$km))),
    pch = shape_list,
    cex = 0.9,
    bty = "n",
    title = expression(K[TC2]^"Cext")
  )
  
  mtext(
    paste0("(", letters[1], ")"),
    side = 3,
    adj = 0.5,
    line = 0.35,
    cex = cex_all,
    font = 2
  )
  

  # experimental data
  #  PTS transporter system (glucose-specific subunit IIA) in Salmonella typhimurium
  par(mar = c(4, 1, 1, 2.5))
  plot(
    mu_rel ~ titrated_level, data = experimental, 
    xlab = NA,
    ylab = NA,
    xlim = xlim, ylim = ylim,
    pch = 19,
    col = uni_lila,
    axes = FALSE,
    cex = cex_all, cex.lab = cex_all
  )
  
  x <- seq(0, 3, length.out = 50)
  y <- 0.15/0.32 + 1.25 * (x / (0.5 + x))^1.2 - 0.22 * x
  lines(x, y, col = uni_lila)
  
  box()
  axis(1)
  axis(2, las=2, labels=FALSE)
  mtext("Titrated level / WT level of transporter", side = 1, outer = FALSE, line = 2.5, cex = cex_all)
  
  legend("bottomright", legend = "Maltose", col = uni_lila, pch = 19, bty = "n", cex=1.05)
  
  
  mtext(
    paste0("(", letters[2], ")"),
    side = 3,
    adj = 0.5,
    line = 0.35,
    cex = cex_all,
    font = 2
  )
  
  dev.off()
}

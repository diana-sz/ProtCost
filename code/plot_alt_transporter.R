library(here)
library(RColorBrewer)

setwd(here())

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



modelnames <- c("A8alt_trans", "A8alt_trans_rev")

for(modelname in modelnames){
  results <- read.csv(paste0("data/", modelname, ".csv"))
  results <- results[results$convergence == 4,]  # plot only simulations that converged
  
  png(paste0("figures/", modelname,".png"),
      type="cairo", units="cm", res = 300,
      width=16, height=7)
  
  par(mar=c(5,4.5,0.5,0.5), mfrow = c(1,2))

  # Set factor levels
  results$kcat <- results$kcat
  results$kcat <- factor(results$kcat, levels = unique(results$kcat))
  results$km <- factor(results$km, levels = unique(results$km))
  opt_phi <- which.max(results[results$kcat == 50.7 & results$km == 1, "mu"])
  results$rel_phi <- results$phi / results$phi[opt_phi]
  results$mu_rel <- results$mu / results$mu[opt_phi]
  
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
    xlab = "Fixed level / optimal level",
    ylab = "Relative growth rate",
    xlim = c(0, 3), ylim = c(0, 1.1),
    pch = point_shapes,
    col = point_colors,
    cex = 1.05, cex.lab = 1.1, cex.axis = 1.05
    #main = "Simulations"
  )
  
  # Add legend
  legend(
    0.7, 0, yjust = 0,
    legend = unique(results$kcat),
    col = palette_colors,
    pch = 19,
    cex = 1.05,
    bty = "n",
    title = expression("k"[ "cat" ])
  )
  
  # Add legend
  legend(
    2, 0, yjust = 0,
    legend = c(unique(results$km)),
    col = rep("black", length(unique(results$km))),
    pch = shape_list,
    cex = 1.05,
    bty = "n",
    title = "Km"
  )
  

  # experimental data
  #  PTS transporter system (glucose-specific subunit IIA) in Salmonella typhimurium
  plot(
    mu_rel ~ titrated_level, data = experimental, 
    xlab = "Titrated level / WT level",
    ylab = "Relative growth rate",
    xlim = c(0, 3), ylim = c(0, 1.1),
    pch = 19,
    col = uni_lila,
    cex = 1.05, cex.lab = 1.1, cex.axis = 1.05
  )
  
  x <- seq(0, 3, length.out = 50)
  y <- 0.15/0.32 + 1.25 * (x / (0.5 + x))^1.2 - 0.22 * x
  lines(x, y, col = uni_lila)
  
  legend("bottomright", legend = "Maltose", col = uni_lila, pch = 19, bty = "n", cex=1.05)
  
  dev.off()
}

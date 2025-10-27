library(here)
library(RColorBrewer)

setwd(here())

modelnames <- c("A8alt_trans", "A8alt_trans_rev")


for(modelname in modelnames){
  results <- read.csv(paste0("data/", modelname, "_alt_trans.csv"))
  results <- results[results$convergence == 4,]  # plot only simulations that converged
  
  png(paste0("figures/", modelname,"_alt_transporter.png"),
      type="cairo", units="cm", res = 300,
      width=9, height=8)
  
  par(mar=c(5,5,2,1))

  # Set factor levels
  results$kcat <- factor(results$kcat, levels = kcats_to_test)
  results$km <- factor(results$km, levels = kms_to_test)
  
  # Set colors based on kcat
  palette_colors <- brewer.pal(length(kcats_to_test), "Paired")
  color_map <- setNames(palette_colors, kcats_to_test)
  point_colors <- color_map[as.character(results$kcat)]
  
  # Set shapes based on km
  shape_list <- 19:(18 - length(kms_to_test) + 1)
  shape_map <- setNames(shape_list, kms_to_test)
  point_shapes <- shape_map[as.character(results$km)]
  
  # Create the plot
  plot(
    mu ~ phi, data = results,
    xlab = "Titrated level / opt. level",
    ylab = "Fraction of optimal growth rate",
    xlim = c(0, 0.2), ylim = c(0, max(results$mu)*1.05),
    pch = point_shapes,
    col = point_colors,
    cex = 1.3, cex.lab = 1.5
  )
  
  # Add legend
  legend(
    "bottomright",
    legend = kcats_to_test,
    col = palette_colors,
    pch = 19,
    cex = 1.3,
    title = "kcat"
  )
  
  # Add legend
  legend(
    "bottomleft",
    legend = c(kms_to_test),
    col = rep("black", length(kms_to_test)),
    pch = shape_list,
    cex = 1.3,
    title = "Km"
  )
  
  dev.off()
}


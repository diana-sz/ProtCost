library(here)
library(RColorBrewer)

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

setwd(here())

cex_all <- 0.95
cex_points <- 0.95
phi_xlim <- 0.4

modelname <- "M10fuel_efflux_rev"

results <- read.csv(here("data", paste0(modelname, "_tradeoff.csv")))
results <- results[results$convergence == 4,]  # plot only simulations that converged
results <- results[results$fuel_phi <= phi_xlim,]
results$efflux_phi <- results$p.EFFLUX/results$c.P

# Set colors based on kcat
kcats <- unique(results$kcat)
#palette_colors <- brewer.pal(length(kcats), "Paired")
palette_colors <- setNames(viridisLite::magma(length(kcats), end = 0.9), kcats)

color_map <- setNames(palette_colors, kcats)

# draws one line per kcat
add_kcat_lines <- function(y_col) {
  for (k in kcats) {
    one_k <- results[results$kcat == k, ]
    one_k <- one_k[order(one_k$fuel_phi), ]
    lines(one_k$fuel_phi, one_k[[y_col]], col = color_map[as.character(k)], lwd = 4)
  }
}

for (format in c("tif", "png")) {
  open_figure("Figure_5", 18, 7.6, format)
  par(mfrow = c(1, 2), oma = c(0, 0, 0.2, 0))
  
  # ---- Panel A: growth rate ----
  par(mar = c(3.7, 3.7, 1, 1))
  
  plot(NA,
       ylab = NA,
       xlab = NA,
       xlim = c(0, phi_xlim),
       ylim = range(results$mu),
       axes = FALSE)
  add_kcat_lines("mu")
  
  box()
  axis(1, cex.axis = cex_all)
  axis(2, las=2, cex.axis = cex_all)
  
  legend("topright", col = rev(palette_colors), 
         legend = rev(kcats), lwd = 4, cex = cex_points, ncol = 2,
         title = expression("efflux pump " * k["cat"]), bty = "n")
  
  mtext(bquote("Growth rate" ~ "[" * h^-1 * "]"), side = 2, outer = FALSE, line = 2.5, cex = cex_all)
  mtext(expression("Proteome fraction fuel reaction " * italic("\u03D5")["FP"]), 
        side = 1, outer = FALSE, line = 2.3, cex = cex_all)
  
  add_panel_label("A", cex = cex_all * 1.2)
  
  # ---- Panel B: optimal phi of efflux pump ----
  par(mar = c(3.7, 3.7, 1, 1))
  
  plot(NA,
       ylab = NA,
       xlab = NA,
       xlim = c(0, phi_xlim),
       ylim = range(results$efflux_phi),
       axes = FALSE)
  add_kcat_lines("efflux_phi")
  
  box()
  axis(1, cex.axis = cex_all)
  axis(2, las=2, cex.axis = cex_all) #, at = seq(0, 0.25, 0.05), labels = c(0.0, NA, 0.1, NA, 0.2, NA))
  
  mtext(expression("Proteome fraction pump " * italic("\u03D5")["EF"]), side = 2, outer = FALSE, line = 2.7, cex = cex_all)
  mtext(expression("Proteome fraction fuel reaction " * italic("\u03D5")["FP"]), side = 1, outer = FALSE, line = 2.3, cex = cex_all)
  
  add_panel_label("B", cex = cex_all * 1.2)
  
  dev.off()
}


# ---- Source data ----
source_dir <- here("figures", "source_data_figures", "Figure 5")

write_source <- function(df, panel, file_name) {
  out_dir <- file.path(source_dir, panel)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(out_dir, file_name), row.names = FALSE)
}

ord <- order(results$kcat, results$fuel_phi)

# Panel A: growth rate vs proteome fraction of the fuel reaction
panel_a <- data.frame(kcat_EF = results$kcat[ord],
                      phi_FP = results$fuel_phi[ord],
                      growth_rate = results$mu[ord])
write_source(panel_a, "5A", "Figure_5A.csv")

# Panel B: proteome fraction of the efflux pump vs proteome fraction of the fuel reaction
panel_b <- data.frame(kcat_EF = results$kcat[ord],
                      phi_FP = results$fuel_phi[ord],
                      phi_EF = results$efflux_phi[ord])
write_source(panel_b, "5B", "Figure_5B.csv")

writeLines(c(
  "Source data for Figure 5",
  "",
  paste0("Simulations: GBA model ", modelname, " (reversible kinetics)."),
  "Only converged simulations with phi_FP <= 0.4 are included, as in the figure.",
  "",
  "5A/Figure_5A.csv",
  "  kcat_EF: turnover number of the efflux pump [1/h] (line color)",
  "  phi_FP: proteome fraction of the fuel-producing reaction FP",
  "  growth_rate: growth rate [1/h]",
  "",
  "5B/Figure_5B.csv",
  "  kcat_EF, phi_FP: as in 5A",
  "  phi_EF: proteome fraction of the efflux pump"
), file.path(source_dir, "README.txt"))
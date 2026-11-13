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

cex_all <- 0.95
cex_points <- 0.8

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

modelname <- "M9_alt_trans_rev"
xlim <- c(0, 3)
ylim <- c(0, 1.1)


for (format in c("tif", "png")) {
  open_figure("Figure_6", 18, 7.6, format)
  par(mfrow = c(1, 2), oma = c(0, 0, 0.2, 0))
  
  # ---- Panel A: simulations ----
  par(mar = c(3.7, 3.7, 1, 1))
  
  results <- read.csv(here("data", paste0(modelname, ".csv")))
  results <- results[results$convergence == 4,]  # plot only simulations that converged
  
  ref <- results[results$kcat == min(results$kcat) & results$km == min(results$km), ]
  ref_row <- ref[which.max(ref$mu), ]
  results$rel_phi <- results$phi / ref_row$phi
  results$mu_rel <- results$mu / ref_row$mu
  results <- results[results$rel_phi <= 3,]
  
  # Set colors based on kcat
  kcat_vals <- sort(unique(results$kcat))
  palette_colors <- viridisLite::magma(length(kcat_vals), end = 0.9)
  color_map <- setNames(palette_colors, kcat_vals)
  
  # Set line types based on km
  km_vals <- sort(unique(results$km))
  lty_list <- seq_along(km_vals)
  lty_map <- setNames(lty_list, km_vals)
  
  plot(NA,
       xlab = NA,
       ylab = NA,
       xlim = xlim,
       ylim = ylim,
       axes = FALSE)
  
  for (kc in kcat_vals) {
    for (km in km_vals) {
      sub <- results[results$kcat == kc & results$km == km, ]
      sub <- sub[order(sub$rel_phi), ]
      lines(sub$rel_phi, sub$mu_rel,
            col = color_map[as.character(kc)],
            lty = lty_map[as.character(km)],
            lwd = 3)
    }
  }
  
  mtext("Relative growth rate", side = 2, outer = FALSE, line = 2.6, cex = cex_all)
  mtext(expression("TS1 proteome fraction relative to optimum " * italic("\u03D5") * "/" * italic("\u03D5")^"\u204E"),
        side = 1, outer = FALSE, line = 2.5, cex = cex_all, adj = 0.7)
  
  box()
  axis(1, cex.axis = cex_all)
  axis(2, las = 2, cex.axis = cex_all)
  
  legend(
    1.2, 0, yjust = 0,
    legend = round(kcat_vals),
    col = palette_colors,
    #lty = 1,
    #lwd = 3,
    pch = 19,
    cex = cex_all,
    bty = "n",
    title = expression(k[cat]^"TS2")
  )
  
  legend(
    2, 0, yjust = 0,
    legend = km_vals,
    col = rep("black", length(km_vals)),
    lty = lty_list,
    lwd = 3,
    cex = cex_all,
    bty = "n",
    title = expression(K[TS2]^"Sext")
  )
  
  add_panel_label("A", cex = cex_all * 1.2)
  
  # ---- Panel B: experimental data ----
  # PTS transporter system (glucose-specific subunit IIA) in Salmonella typhimurium
  par(mar = c(3.7, 3.7, 1, 1))
  
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
  axis(2, las = 2, cex.axis = cex_all)
  mtext("Relative growth rate", side = 2, outer = FALSE, line = 2.6, cex = cex_all)
  mtext("Titrated level / WT level of transporter", side = 1, outer = FALSE, line = 2.5, cex = cex_all)
  
  legend("bottomright", legend = "Maltose", col = "black", pch = 19, bty = "n", cex = cex_all)
  
  add_panel_label("B", cex = cex_all * 1.2)
  
  dev.off()
}


# ---- Source data ----
source_dir <- here("figures", "source_data_figures", "Figure 6")

write_source <- function(df, panel, file_name) {
  out_dir <- file.path(source_dir, panel)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(out_dir, file_name), row.names = FALSE)
}

# Panel A: simulations (color = kcat of TS2, line type = Km of TS2)
panel_a <- data.frame(kcat_TS2 = results$kcat,
                      Km_TS2 = results$km,
                      phi_TS1_relative = results$rel_phi,
                      relative_growth_rate = results$mu_rel)
panel_a <- panel_a[order(panel_a$kcat_TS2, panel_a$Km_TS2,
                         panel_a$phi_TS1_relative), ]
write_source(panel_a, "6A", "Figure_6A.csv")

# Panel B: experimental points and phenomenological fit
panel_b_exp <- data.frame(titrated_level = experimental$titrated_level,
                          relative_growth_rate = experimental$mu_rel)
write_source(panel_b_exp, "6B", "Figure_6B_van_der_Vlag_1994.csv")

panel_b_fit <- data.frame(titrated_level = x,
                          relative_growth_rate = y)
write_source(panel_b_fit, "6B", "Figure_6B_fit_Bruggeman_2023.csv")

writeLines(c(
  "Source data for Figure 6",
  "",
  paste0("6A/Figure_6A.csv: simulations with GBA model ", modelname,
         " (reversible kinetics)."),
  "  Only converged simulations with phi_TS1_relative <= 3 are included, as in the figure.",
  "  kcat_TS2: turnover number of the alternative transporter TS2 [1/h] (line color)",
  "  Km_TS2: Michaelis constant of TS2 for S_ext [g/L] (line type)",
  "  phi_TS1_relative: proteome fraction of TS1 relative to its optimum (phi/phi*)",
  "  relative_growth_rate: growth rate relative to the optimum",
  "",
  "6B/Figure_6B_van_der_Vlag_1994.csv: experimental data (points) from",
  "  van der Vlag et al. (1994), redrawn from Bruggeman et al. (2023).",
  "  PTS IIA(Glc) in Salmonella typhimurium grown on maltose.",
  "  titrated_level: titrated level / wild-type level of the transporter",
  "  relative_growth_rate: growth rate relative to the wild type",
  "",
  "6B/Figure_6B_fit_Bruggeman_2023.csv: phenomenological fit (line) from",
  "  Bruggeman et al. (2023):",
  "  relative_growth_rate = 0.15/0.32 + 1.25 * (x / (0.5 + x))^1.2 - 0.22 * x,",
  "  with x = titrated_level."
), file.path(source_dir, "README.txt"))
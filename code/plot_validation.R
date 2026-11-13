library(RColorBrewer)
library(scales)
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

filelist <- c("M8_GBA.csv",
              "M8_rev_GBA.csv")

cex_all <- 1.1
colors <- brewer.pal(8, "Paired")[c(2,6)]

plot_and_fit <- function(data, growth_rates, ylim, xlim,
                         show_yaxis = FALSE,
                         xlab = FALSE,
                         ylab = FALSE,
                         legend = TRUE,
                         fit = TRUE) {
  
  if (length(data) > length(colors)) {
    print("Too many groups")
    return()
  }
  
  for (i in seq_along(data)) {
    df <- data.frame(x = growth_rates, y = data[[i]])
    
    if (i == 1) {
      plot(
        y ~ x, data = df,
        col = colors[i],
        pch = 20,
        ylim = ylim,
        xlim = xlim,
        axes = FALSE,
        xlab = "",
        ylab = "",
        cex = cex_all
      )
    } else {
      points(df$x, df$y, col = colors[i], pch = 20, cex = cex_all)
    }
    
    if (fit) {
      fitted <- lm(y ~ x, data = df)
      x_seq <- seq(min(xlim), max(xlim), length.out = 200)
      pred <- predict(fitted,
                      newdata = data.frame(x = x_seq),
                      interval = "confidence")
      
      polygon(
        c(x_seq, rev(x_seq)),
        c(pred[, "lwr"], rev(pred[, "upr"])),
        col = alpha(colors[i], 0.3),
        border = NA
      )
      lines(x_seq, pred[, "fit"], col = colors[i], lwd = 1.5)
    }
  }
  
  if (legend) {
    legend(
      "topleft",
      legend = names(data),
      col = colors[seq_along(data)],
      pch = 20,
      bty = "n",
      cex = cex_all*0.9
    )
  }
}

for (format in c("tif", "png")) {
  open_figure("Appendix_Figure_S4", 18, 8, format)
  par(mfrow=c(1,2), oma = c(1,1.5,0.5,0))


  for (n in seq_along(filelist)) {
    left <- ifelse (n == 1, 2.5, 1)
    par(mar = c(2.5, left, 1, 3.5-left))
  
    filename <- filelist[n]
    modelname <- strsplit(filename, "_GBA.csv")[[1]][1]
  
    opt_data <- read.csv(here("data", filename), row.names = 1)
  
    reactants <- colnames(opt_data)[5:which(colnames(opt_data) == "P")]
    concentrations <- opt_data[, reactants]
    p_opt <- opt_data[, grep("p.", colnames(opt_data))]
    phi_opt <- p_opt / opt_data[, "P"]
  
    rnas <- rowSums(concentrations[, grep("RNA", colnames(concentrations)), drop = FALSE])
  
    sectors <- list(
      "Translation" = phi_opt$p.R,
      "RNA/protein ratio" = rnas / opt_data[, "P"]
    )
  
    plot_and_fit(
      sectors,
      opt_data$mu,
      ylim = c(0, 0.75),
      xlim = c(0, max(opt_data$mu)),
      show_yaxis = (n == 1),
      legend = (n == 1)
    )
  
    # Panel title (centered)
    label <- ifelse(grepl("rev", filename), "reversible", "irreversible")
    mtext(
      paste(label, "kinetics"),
      side = 3,
      adj = 0.5,
      line = 0.35,
      cex = cex_all*0.9,
      font = 2
    )
  
    # Panel letter
    add_panel_label(LETTERS[n], cex = 0.97 * 1.2)
  
    labels <- c("0.0", "", "0.2", "", "0.4", "", "0.6", "")
    if (n == 2) labels <- FALSE
    axis(1)
    axis(2, las = 2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), labels = labels)
    box()
  
  }

  mtext(bquote("Growth rate" ~ "[" * h^-1 * "]"),
        side = 1, outer = TRUE, line = 0, cex = cex_all)
  mtext("Proteome fraction / ratio",
        side = 2, outer = TRUE, line = 0.05, cex = cex_all)

  dev.off()
}

library(RColorBrewer)
library(scales)
library(here)

setwd(paste0(here(), "/data"))

filelist <- c("M8_GBA.csv",
              "M8_rev_GBA.csv")

cex_all <- 1.1
colors <- brewer.pal(8, "Dark2")

plot_and_fit <- function(data, growth_rates, ylim, xlim,
                         show_yaxis = FALSE,
                         xlab = FALSE,
                         ylab = FALSE,
                         legend = TRUE,
                         fit = TRUE) {
  
  if (length(data) > 8) {
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

png(
  paste0("../figures/M8_RP_ratio.png"),
  type = "cairo",
  units = "cm",
  width = 18,
  height = 8,
  res = 300
)

par(mfrow=c(1,2), oma = c(1.5,1.5,0.5,0))


for (n in seq_along(filelist)) {
  left <- ifelse (n == 1, 2.5, 1)
  par(mar = c(2.5, left, 1, 3.5-left))
  
  filename <- filelist[n]
  modelname <- strsplit(filename, "_GBA.csv")[[1]][1]
  
  opt_data <- read.csv(filename, row.names = 1)
  
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
    ylim = c(0, 0.6),
    xlim = c(0, max(opt_data$mu)),
    show_yaxis = (n == 1),
    legend = (n == 1)
  )
  
  label <- ifelse(grepl("rev", filename), "reversible", "irreversible")
  mtext(
    paste0("(", letters[n], ") ", label),
    side = 3,
    adj = 0.5,
    line = 0.35,
    cex = cex_all,
    font = 2
  )
  
  labels <- c("0.0", "", "0.2", "", "0.4", "", "0.6")
  if (n == 2) labels <- FALSE
  axis(1)
  axis(2, las = 2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = labels)
  box()
  
}

mtext(bquote("Growth rate"~ mu ~ (h^-1)),
      side = 1, outer = TRUE, line = 0, cex = cex_all)
mtext("Proteome fraction / ratio",
      side = 2, outer = TRUE, line = -0.2, cex = cex_all)

dev.off()

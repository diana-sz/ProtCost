library(RColorBrewer)
library(scales)
library(here)

filelist <- c("M8_GBA.csv",
              "M8_rev_GBA.csv")

cex_all <- 1.05
colors <- brewer.pal(8, "Paired")[c(2,6,7)]

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

png(
  here("figures", "M8_prot_laws.png"),
  type = "cairo",
  units = "cm",
  width = 19,
  height = 8,
  res = 300
)

par(mfrow=c(2,4), oma = c(2.5,2.5,0,0))

sectors_to_plot <- list(
  "1" = c(
    TS = "p.TS",
    AAS = "p.AAS"
  ),
  "2" = c(
    ATPS = "p.ATPS",
    NTS = "p.NTS"
  ),
  "3" = c(
    ADPS = "p.ADPS",
    DNAP = "p.DNAP"
  ),
  "4" = c(
    R = "p.R",
    RNAP = "p.RNAP"
  )
)


for (n in seq_along(filelist)) {
  par(mar = c(1.3, 1.2, 2.2, 2))
  
  filename <- filelist[n]
  modelname <- strsplit(filename, "_GBA.csv")[[1]][1]
  
  opt_data <- read.csv(here("data", filename), row.names = 1)
  
  reactants <- colnames(opt_data)[5:which(colnames(opt_data) == "P")]
  concentrations <- opt_data[, reactants]
  p_opt <- opt_data[, grep("p.", colnames(opt_data))]
  phi_opt <- p_opt / opt_data[, "P"]
  
  label <- ifelse(grepl("rev", filename), "reversible", "irreversible")
  
  for (s in names(sectors_to_plot)){
    sector_names <- sectors_to_plot[[s]]
    sectors <- list()
    ylim <- 0
    for (name in names(sector_names)){
      sectors[name] = phi_opt[sector_names[[name]]]
      ylim <- max(c(max(phi_opt[sector_names[[name]]]), ylim))
    }
    
    
    plot_and_fit(
      sectors,
      opt_data$mu,
      ylim = c(0, ylim*1.1),
      xlim = c(0, 1.1),
      show_yaxis = TRUE,
      legend = (n == 1)
    )

  axis(1, cex.axis = 0.85, padj=-1)
  axis(2, las = 2, cex.axis = 0.85, hadj=0.85)
  box()
  
  }
  
  mtext(
    label,
    side = 3,
    outer = FALSE,
    line = 0.5,
    cex = cex_all,
    at = -1.75
  )
}


mtext(bquote("Growth rate" ~ "[" * h^-1 * "]"),
      side = 1, outer = TRUE, line = 1.3, cex = cex_all)
mtext(expression("Proteome fraction " * italic("\u03A6")), 
      side = 2, outer = TRUE, line = 1, cex = cex_all)


dev.off()

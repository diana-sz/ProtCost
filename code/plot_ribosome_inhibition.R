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

cex_all <- 1.1
modelnames <- c("M8_inh",
                "M8_inh_rev")

for (format in c("tif", "png")) {
  open_figure("Appendix_Figure_S5", 18, 8, format)
  par(mfrow=c(1,2), oma = c(1,1.5,0.5,0))

  for (n in seq_along(modelnames)) {
    left <- ifelse (n == 1, 2.5, 1)
    par(mar = c(2.5, left, 1, 3.5-left))
  
    modelname <- modelnames[n]
    opt_data <- read.csv(here("data", paste0(modelname, ".csv")), row.names = 1)
  
    # filter out results that did not converge / do not grow
    opt_data <- opt_data[opt_data$convergence == 4, ]
    opt_data <- opt_data[opt_data$mu > 1e-6, ]
  
    opt_data$phiR <- opt_data$p.R/rowSums(opt_data[, grep("p\\.", colnames(opt_data))])
    palette_colors <- viridisLite::magma(length(unique(opt_data$x_C)), end = 0.9)
    color_map <- setNames(palette_colors, unique(opt_data$x_C))
    opt_data$color <- color_map[as.character(opt_data$x_C)]
  
    for(x_c in unique(opt_data$x_C)){
      one_w <- opt_data[opt_data$x_C == x_c, ]
    
      if(nrow(one_w) < 3) next
    
      plot(phiR ~ mu, data=one_w, col = one_w$color,
           xlim = c(0,1.1), ylim = c(0, 0.9),
           xlab = NA, ylab = NA,
           axes = FALSE,
           pch = 19)
    
      fit <- lm(phiR ~ mu, data = one_w)
      mu_range <- range(0, max(one_w$mu, na.rm = TRUE))
      mu_seq <- seq(mu_range[1], mu_range[2], length.out = 100)
      rp_pred <- predict(fit, newdata = data.frame(mu = mu_seq))
      lines(mu_seq, rp_pred, col = unique(one_w$color), lty = 1)
    
      par(new=TRUE)
    }
    par(new=FALSE)
  
    # Panel title (centered)
    label <- ifelse(grepl("rev", modelname), "reversible", "irreversible")
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
  
    labels <- ifelse(n == 2, FALSE, TRUE)
    axis(1)
    axis(2, las = 2, labels = labels)
    box()
  
    if(n == 1){
      mtext(expression("Ribosomal proteome fraction " * italic("\u03D5")["R"]), 
            side = 2, outer = FALSE, line = 2.3, cex = cex_all*0.9)
    }
  
    no_inh <- opt_data[opt_data$x_INH== min(opt_data$x_INH),]
    fit_no_inh <- lm(phiR ~ mu, data = no_inh)
    abline(fit_no_inh, lty = 2)
  
    ### arrows ####
    arrow_length <- 0.3
  
    # arrow for glucose concentration
    slope <- coef(fit_no_inh)[2]
  
    # Choose a point above the fitted line to start the arrow
    x_end <- 0.8
    y_end <- predict(fit_no_inh, newdata = data.frame(mu = x_end)) - 0.05
  
    # Define arrow length
    x_start <- x_end - arrow_length
    y_start <- y_end - arrow_length * slope
  
    # Add arrow and label
    arrows(x0 = x_start, y0 = y_start, x1 = x_end, y1 = y_end,
           length = 0.1, col = "black", lwd = 1.5)
  
    text(x = 0.68, y = 0.2, labels = "Substrate\nconcentration",
         adj = c(0, 0.5), font = 1, cex = cex_all*0.8)
  
  
    # arrow for inhibitor concentration
    max_xc <- max(opt_data$x_C)
    highlight_group <- opt_data[opt_data$x_C == max_xc, ]
  
    # Fit linear model
    fit_high <- lm(phiR ~ mu, data = highlight_group)
    slope <- coef(fit_high)[2]
  
    # Choose a point above the fitted line to start the arrow
    x_start <- 0.9
    y_start <- predict(fit_high, newdata = data.frame(mu = x_start)) + 0.05
  
    # Define arrow length
    x_end <- x_start - arrow_length
    y_end <- y_start - arrow_length * slope
  
    # Add arrow and label
    arrows(x0 = x_start, y0 = y_start, x1 = x_end, y1 = y_end,
           length = 0.1, col = "black", lwd = 1.5)
  
    text(x = 0.68, y = 0.8, labels = "R inhibitor\nconcentration",
         adj = c(0, 0.5), font = 1, cex = cex_all*0.8)
  }

  mtext(bquote("Growth rate" ~ "[" * h^-1 * "]"),
        side = 1, outer = TRUE, line = 0, cex = cex_all)
  dev.off()
}
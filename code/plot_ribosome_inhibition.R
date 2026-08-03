library(RColorBrewer)
library(scales)
library(here)

cex_all <- 1.1
modelnames <- c("M8_inh",
                "M8_inh_rev")



png(here("figures", "M8_inh_ribosome_inhibition.png"), 
    type="cairo", units="cm", res = 300,
    width = 18, height = 8)

par(mfrow=c(1,2), oma = c(1.5,1.5,0.5,0))

for (n in seq_along(modelnames)) {
  left <- ifelse (n == 1, 2.5, 1)
  par(mar = c(2.5, left, 1, 3.5-left))
  
  modelname <- modelnames[n]
  opt_data <- read.csv(here("data", paste0(modelname, ".csv")), row.names = 1)
  opt_data <- opt_data[opt_data$convergence == 4, ]
  
  opt_data$phiR <- opt_data$p.R/rowSums(opt_data[, grep("p\\.", colnames(opt_data))])
  palette_colors <- brewer.pal(length(unique(opt_data$x_C)), "Paired")
  color_map <- setNames(palette_colors, unique(opt_data$x_C))
  opt_data$color <- color_map[as.character(opt_data$x_C)]
  
  for(x_c in unique(opt_data$x_C)){
    one_w <- opt_data[opt_data$x_C == x_c, ]
    
    if(all(one_w$mu < 1e-6)) next
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
  
  label <- ifelse(grepl("rev", modelname), "reversible", "irreversible")
  mtext(
    paste0("(", letters[n], ") ", label),
    side = 3,
    adj = 0.5,
    line = 0.35,
    cex = cex_all,
    font = 2
  )
  
  labels <- ifelse(n == 2, FALSE, TRUE)
  axis(1)
  axis(2, las = 2, labels = labels)
  box()
  
  if(n == 1){
    mtext(expression("Ribosomal proteome fraction " * italic("\u03A6")["R"]), 
          side = 2, outer = FALSE, line = 2.6, cex = cex_all)
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
  
  text(x = 0.68, y = 0.8, labels = "R inhibitior\nconcentration",
       adj = c(0, 0.5), font = 1, cex = cex_all*0.8)
}

mtext(bquote("Growth rate" ~ "[" * h^-1 * "]"),
      side = 1, outer = TRUE, line = 0, cex = cex_all)
dev.off()


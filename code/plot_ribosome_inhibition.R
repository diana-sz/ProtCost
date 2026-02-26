library(RColorBrewer)
library(scales)
library(here)

cex_all <- 1

modelname <- "M8_inh"
opt_data <- read.csv(here("data", paste0(modelname, ".csv")), row.names = 1)
opt_data <- opt_data[opt_data$convergence == 4, ]

opt_data$phiR <- opt_data$p.R/rowSums(opt_data[, grep("p\\.", colnames(opt_data))])
palette_colors <- brewer.pal(length(unique(opt_data$x_C)), "Paired")
color_map <- setNames(palette_colors, unique(opt_data$x_C))
opt_data$color <- color_map[as.character(opt_data$x_C)]


png(here("figures", paste0(modelname, "_ribosome_inhibition.png")), 
    type="cairo", units="cm", res = 300,
    width=9, height=7.5)

par(mfrow=c(1,1), oma = c(0,0,0.5,0))
par(mar = c(3.7, 3.7, 0.5, 1))

for(x_c in unique(opt_data$x_C)){
  one_w <- opt_data[opt_data$x_C == x_c, ]
  if(nrow(one_w) < 3){break}
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

axis(1)
axis(2, las = 2)
box()


mtext(expression("Ribosomal proteome fraction " * Phi["R"]), 
      side = 2, outer = FALSE, line = 2.6, cex = cex_all)
mtext(bquote("Growth rate" ~ "[" * h^-1 * "]"), side = 1, outer = FALSE, line = 2.5, cex = cex_all)


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

text(x = x_end - 0.1, y = y_end - 0.15, labels = "Substrate\nconcentration",
     adj = c(0, 0.5), font = 1, cex = cex_all*0.9)


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

text(x = x_end + 0.05, y = y_end + 0.1, labels = "R inhibitior\nconcentration",
     adj = c(0, 0.5), font = 1, cex = cex_all*.9)

dev.off()


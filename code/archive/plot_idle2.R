library(here)
library(RColorBrewer)

xlim <- c(0, 0.45)
ylim <- c(0, 0.45)
xlab_line <- 2.15
cex_all <- 0.97

for(modelname in c("M10_IDLE2")){
  data <- read.csv(here("data", paste0(modelname, ".csv")))
  data <- data[data$convergence == 4, ]

  png(here("figures", paste0(modelname, ".png")), 
      type="cairo", units="cm",
      width=9, height=7.5, res=300)
  par(mfcol=c(1,1), mar = c(3,3.5,0.1,0.1))
  
  mu <- round(data$mu, 2)
  cols <- setNames(rev(brewer.pal(length(unique(mu)), "Paired")), unique(mu))
  point_cols <- cols[as.character(mu)] 

  plot(data$idle1_phi, data$idle2_phi,
       xlim = xlim, ylim = ylim,
       col = point_cols,
       pch=20,
       axes = FALSE,
       xlab = NA,
       ylab=NA)
  
  # linear fits per mu
  for (m in unique(mu)) {
    idx <- mu == m
    fit <- lm(idle1_phi ~ idle2_phi, data = data[idx, ])
    abline(fit, col = cols[as.character(m)], lwd = 1.5)
  }
  
  box()
  axis(1, cex.axis=cex_all)
  axis(2,las=2, cex.axis=cex_all)
  mtext(bquote("IDLE1 proteome fraction " * italic("\u03A6")["IDLE1"]), side = 1, cex = cex_all, line = xlab_line, outer = FALSE)
  mtext(bquote("IDLE2 proteome fraction " * italic("\u03A6")["IDLE2"]), side = 2, cex = cex_all, line = xlab_line+0.1, outer = FALSE)
  
  legend("topright", legend = unique(mu), col =  cols, pch = 20,
         title = "Growth rate", cex = 0.9)
    
  dev.off()
}


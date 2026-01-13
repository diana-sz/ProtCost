library(here)
library(RColorBrewer)

setwd(here())

xlim <- c(0, 0.45)
ylim <- c(0, 0.45)
xlab_line <- 2
cex_all <- 1.05

for(modelname in c("M10_GFP_RFP_rev")){
  data <- read.csv(paste0("data/", modelname, ".csv"), row.names = 1)
  data <- data[data$convergence == 4, ]
  data <- data[complete.cases(data),]
  
  png(paste0("figures/", modelname, ".png"), 
      type="cairo", units="cm",
      width=9, height=7, res=300)
  par(mfcol=c(1,1), mar = c(3,3.5,0.1,0.1))
  
  mu <- round(data$mu, 2)
  cols <- setNames(rev(brewer.pal(length(unique(mu)), "Paired")), unique(mu))
  point_cols <- cols[as.character(mu)] 

  plot(data$gfp_phi, data$rfp_phi,
       xlim = xlim, ylim = ylim,
       col = point_cols,
       pch=20,
       axes = FALSE,
       xlab = NA,
       ylab=NA)
  
  # linear fits per mu
  for (m in unique(mu)) {
    idx <- mu == m
    fit <- lm(rfp_phi ~ gfp_phi, data = data[idx, ])
    abline(fit, col = cols[as.character(m)], lwd = 1.5)
  }
  
  box()
  axis(1)
  axis(2,las=2)
  mtext(bquote("GFP proteome fraction"), side = 1, cex = cex_all, line = xlab_line, outer = FALSE)
  mtext(bquote("RFP proteome fraction"), side = 2, cex = cex_all, line = xlab_line+0.3, outer = FALSE)
  
  legend("topright", legend = unique(mu), col =  cols, pch = 20,
         title = "Growth rate")
    
  dev.off()
}





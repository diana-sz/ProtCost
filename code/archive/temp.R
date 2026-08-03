library(here)
library(RColorBrewer)

cex_all <- 0.97
cex_points <- 0.95

xlim <- c(0, 3)
ylim <- c(0, 1.1)

modelname <- "M8"
results <- read.csv(here("data", paste0(modelname, "_protein_cost.csv")))
results <- results[results$convergence == 4,]  # plot only simulations that converged
results <- results[results$reaction == "TS", ]  

png(here("figures", "TS1.png"),
    type="cairo", units="cm", res = 300,
    width=9, height=7.5)


par(mfrow=c(1,1), oma = c(0,0,0,0))
par(mar = c(3.7, 3.7, 0.5, 1))

plot(
  mu_norm ~ rel_phi, data = results,
  xlab = NA,
  ylab = NA,
  xlim = xlim,
  ylim = ylim,
  pch = 19,
  col = "black",
  axes = FALSE
)


box()
axis(1, cex.axis = cex_all)
axis(2, las=2, cex.axis = cex_all)
mtext("Relative growth rate", side = 2, outer = FALSE, line = 2.6, cex = cex_all)
mtext("Relative proteome fraction", side = 1, outer = FALSE, line = 2.5, cex = cex_all)



dev.off()


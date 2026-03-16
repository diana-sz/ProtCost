library(here)
library(RColorBrewer)

modelname <- "M8_rev"
tested_kcat_factors <- c("1", "2", "5", "100", "10000")

all_data <- data.frame()

for(factor in tested_kcat_factors){
  modelname <- ifelse(factor == "10000", "M8", "M8_rev")
  data <- read.csv(here("data", paste0(modelname, "_prot_cost_kf_", factor,".csv")))
  data$factor <- factor
  all_data <- rbind(all_data, data)
  
  print( max(data$mu_norm))

}

colors <- brewer.pal(length(tested_kcat_factors), "Paired")
names(colors) <- tested_kcat_factors

for(reaction in unique(all_data$reaction)){
  one_prot <- all_data[all_data$reaction == reaction,]

  plot(NA, xlim = c(0,max(one_prot$phi)), ylim = c(0,1), xlab = NA, ylab = NA)
  
  points(one_prot$phi, one_prot$mu, col = colors[one_prot$factor], pch = 19)
  
  text(0.4, 0.25, reaction, cex=1.4)
  
  tested_kcat_factors[length(tested_kcat_factors)] <- "irreversible"
  legend("bottomright", tested_kcat_factors, pch = 19, col = colors)
  
  mtext(expression("Growth rate [h" ^ -1 * "]" ), side=2, cex = 1.2, line = 2.4)
  mtext(expression("Proteome fraction " * italic("\u03A6")), side=1, cex = 1.2, line = 2.4)
  
}


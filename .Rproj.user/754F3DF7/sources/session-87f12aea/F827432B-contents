library(RColorBrewer)
library(scales)
library(here)

setwd(paste0(here(), "/data"))

filelist <- c("A7simple_GBA.csv",
              "A7simple_rev_GBA.csv")

size <- 1.5

plot_and_fit <- function(data, growth_rates, ylim, xlim, 
                         xlab=bquote("Growth rate"~ mu ~ (h^-1)),
                         ylab="Proteome fraction / ratio", legend=TRUE,fit=TRUE){
  if(length(data) > 8){
    print("Too many groups")
    return()
  }
  
  colors <- brewer.pal(8, "Dark2")
  for(i in seq_along(data)){
    df <- data.frame(x=growth_rates, y=data[[names(data)[i]]])
    plot(y~x, data=df,
         col = colors[i],
         pch=20,
         ylim=ylim, 
         xlim=xlim,
         ylab=ylab,
         xlab=xlab,
         cex=1.5,
         cex.lab = 1.3)
    
    if(fit){
      fitted <- lm(y~x, data=df)
      
      x_seq <- seq(-0.1,1.3,0.01)
      pred <- predict(fitted, newdata = data.frame(x = x_seq), interval = "confidence")
      
      # Add shaded confidence interval
      polygon(c(x_seq, rev(x_seq)), 
              c(pred[, "lwr"], rev(pred[, "upr"])), 
              col = alpha(colors[i],0.3), border = NA)  # Transparent shading
      abline(fitted, col = colors[i])
    }
    par(new=TRUE)
    
  }
  if(legend){
    legend("topleft", 
           bty = "n",
           legend = names(data),
           col = colors[1:length(data)], pch=20, cex=1.1)
  }
  par(new=FALSE)
}


for (filename in filelist){
  modelname <- strsplit(filename, "_GBA.csv")[[1]][1]
  opt_data <- read.csv(filename, row.names = 1)
  reactants <- colnames(opt_data)[5:which(colnames(opt_data)=="P")]
  concentrations <- opt_data[, reactants]
  p_opt <- opt_data[,grep("p.",colnames(opt_data))]
  phi_opt <- p_opt/opt_data[,"P"]
  ordered_mu_ind <- order(opt_data$mu)
  ordered_mu <- opt_data$mu[ordered_mu_ind]
  # phi_opt_ordered <- phi_opt[ordered_mu_ind,]
  
  #### Make plots ################################################################
  legend <- TRUE

  rnas <- rowSums(concentrations[, grep("RNA", colnames(concentrations)), drop = FALSE])
  sectors <- list("Translation" = phi_opt$p.R,
                  "RNA/protein ratio" = rnas/opt_data[,"P"])
  png(paste0("../figures/", modelname, "_RP_ratio.png"),
      type="cairo", units="cm", 
      width=10, height=9, res=300)
  par(mar=c(4.5,4.5,1,1))
  plot_and_fit(sectors, 
               opt_data$mu, 
               ylim = c(0,0.6), 
               xlim = c(0, max(ordered_mu)),
               legend = legend)
  #title("Predicted")
  dev.off()
  
  
  # biomass composition
  concentrations_ordered <- concentrations[ordered_mu_ind,]
  prot <- ncol(concentrations)
  RNA <- c(grep("RNA", reactants), grep("TC", reactants))
  DNA <- grep("DNA", reactants)
  
  prot_sum <- concentrations_ordered[, prot, drop = FALSE]
  RNA_sum <- rowSums(concentrations_ordered[, unique(RNA), drop = FALSE])
  DNA_sum <- concentrations_ordered[, DNA, drop = FALSE]
  rest_sum <- rowSums(concentrations_ordered[, -c(DNA, RNA, prot), drop = FALSE])
  
  sectors <- cbind(rest_sum, DNA_sum, RNA_sum, prot_sum)
  sectors <-  t(apply(sectors[,ncol(sectors):1], 1, cumsum))
  sector_names <- c("Rest", "DNA", "RNA", "Protein")

  colors <- brewer.pal(ncol(sectors), "Paired")
  rel_comp <- sectors/sectors[, ncol(sectors)]
  
  png(paste0("../figures/", modelname, "_biomass.png"), 
      type="cairo", units="cm",
      width=15, height=12, res=300)
  par(mar = c(5,5,1,1), mfrow=c(1,1))
  plot(NA,
       xlim=c(0, max(ordered_mu)),
       ylim=c(0, 1),#max(rel_comp)),
       yaxs="i", xaxs="i",
       # xaxt = "n", yaxt = "n",
       cex.lab = size, cex.axis = size,
       xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
       ylab="Predicted biomass")
  for(i in ncol(rel_comp):1){
    col <- colors[i]
    polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(rel_comp)), rev(rel_comp[,i])),
            col=col, border=col)
    par(new=TRUE)
  }
  par(new=FALSE)
  legend('bottomleft', sector_names, col = rev(colors), pch = 15, cex=0.8, bg="white")
  dev.off()
  
}




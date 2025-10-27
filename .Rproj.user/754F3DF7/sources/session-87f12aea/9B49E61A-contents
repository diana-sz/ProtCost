library(RColorBrewer)
library(scales)
library(here)

setwd(here())

source("code/process_ecocyc_annotation2.R")
source("code/process_proteomics.R")

scaled <- FALSE

all_row_names <- unique(c(rownames(mori_proteomics), rownames(mori_proteomics2), rownames(schmidt_proteomics)))

# Align each data frame to the combined row names, filling missing rows with NA
df1_aligned <- data.frame(row.names = all_row_names, mori_proteomics[all_row_names, , drop = FALSE])
df2_aligned <- data.frame(row.names = all_row_names, mori_proteomics2[all_row_names, , drop = FALSE])
df3_aligned <- data.frame(row.names = all_row_names, schmidt_proteomics[all_row_names, , drop = FALSE])

# Combine the data frames
merged <- cbind(df1_aligned, df2_aligned, df3_aligned)

proteomics_datasets <- list(Mori_titrated_glucose=mori_proteomics,
                            Mori_c_sources=mori_proteomics2,
                            Schmidt=schmidt_proteomics,
                            Merged=merged)


if(scaled){
  scaled_1 <- mori_growth_rates/max(mori_growth_rates)
  scaled_2 <- mori_growth_rates2/max(mori_growth_rates2)
  scaled_3 <- schmidt_growth_rates/max(schmidt_growth_rates)
  all_growth_rates <- c(scaled_1, scaled_2, scaled_3)
  growth_datasets <- list(Mori_titrated_glucose=scaled_1,
                          Mori_c_sources=scaled_2,
                          Schmidt=scaled_3,
                          Merged=all_growth_rates)
}else{
  all_growth_rates <- c(mori_growth_rates, mori_growth_rates2, schmidt_growth_rates)
  growth_datasets <- list(Mori_titrated_glucose=mori_growth_rates,
                          Mori_c_sources=mori_growth_rates2,
                          Schmidt=schmidt_growth_rates,
                          Merged=all_growth_rates)
}

#### Function definitions ######################################################
get_sector_fraction <- function(targets, proteomics, gene_annot){
  sector_sum <- data.frame()
  
  for (target in targets){
    target_genes <- c()
    for (gene in names(gene_annot)){
      if (target %in% gene_annot[[gene]]){
        target_genes <- c(target_genes, gene)
        sharing_factor <- 1/length(gene_annot[[gene]]) # if gene belongs multiple groups
      }
    }
    if(is.null(target_genes)){
      print(paste("no genes in", target, "group"))
      next
    }
    target_rows <- proteomics[target_genes, ]*sharing_factor
    target_sums <- colSums(target_rows, na.rm = TRUE)
    sector_sum <- rbind(sector_sum, target_sums) 
  }
  
  summed_groups <- colSums(sector_sum, na.rm=TRUE)
  fractions <- summed_groups/colSums(proteomics, na.rm=TRUE)
  names(fractions) <- colnames(proteomics)
  return(fractions)
}


get_sector_fraction_gene_list <- function(target_genes, proteomics){
  target_rows <- proteomics[target_genes, ]
  target_sums <- colSums(target_rows, na.rm = TRUE)
  fractions <- target_sums/colSums(proteomics, na.rm=TRUE)
  names(fractions) <- colnames(proteomics)
  return(fractions)
}


plot_and_fit <- function(data, growth_rates, ylim, xlim, 
                         xlab=bquote("Growth rate"~ mu ~ (h^-1)),
                         ylab="Proteome fraction", legend=TRUE){
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
    
    fit <- lm(y~x, data=df)
    
    x_seq <- seq(-0.1,1.3,0.01)
    pred <- predict(fit, newdata = data.frame(x = x_seq), interval = "confidence")
    
    if(length(growth_rates)>30){
      print(paste0(names(data)[i], " phi at mu=1: ", round(pred[which(x_seq==1)],4),
                  ", mu=0: ", round(pred[which(x_seq==0)],4)))
    }
  
    # Add shaded confidence interval
    polygon(c(x_seq, rev(x_seq)), 
            c(pred[, "lwr"], rev(pred[, "upr"])), 
            col = alpha(colors[i],0.3), border = NA)  # Transparent shading
    abline(fit, col = colors[i])
    
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


#### Make plots ################################################################
upper_groups <- c("glycolytic process",
                  "D-glucose transmembrane transport",
                  "ketone body metabolic process",
                  "energy derivation by oxidation of organic compounds")
#lower_groups <- c("electron transport chain")
met_groups <- c("lipid biosynthetic process",            
                "amino acid biosynthetic process",
                "nucleobase biosynthetic process",
                "carboxylic acid transport",
                "NADPH regeneration",
                upper_groups) 
info_group <- c("DNA-templated transcription",
                "DNA biosynthetic process",
                "RNA processing", "translation",
                "protein maturation",
                "rRNA catabolic process",
                "mRNA catabolic process",
                "tRNA decay")
xlim <- c(0,1.15)

#+ fig.height = 29.7, fig.width = 12

pdf("figures/proteomics_ecocyc.pdf")
all_rib_data <- data.frame()
par(mfrow=c(2, 2), mar = c(5,5,3,1))
#par(mfcol=c(2,2))

for(dataset in names(proteomics_datasets)){
  legend <- TRUE
  #if(dataset != "Merged"){legend <- FALSE}
  current_dataset <- proteomics_datasets[[dataset]]
  current_growth <- growth_datasets[[dataset]]
  
  # big sectors
  Metabolism <- get_sector_fraction(met_groups, current_dataset, gene_annot)
  Information_processing <- get_sector_fraction(info_group, current_dataset, gene_annot)
  sectors <- list("Metabolism" = Metabolism,
                  "Genetic Information Processing" = Information_processing)
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.5), 
               xlim = xlim,
               legend = legend)
  title(dataset)
  
  sectors <- data.frame(cbind(Information_processing, Metabolism))
  sectors$Other <- 1-rowSums(sectors)
  ordered_mu <- current_growth[order(current_growth)]
  ordered_phis <- sectors[order(current_growth),]
  colors <- brewer.pal(3, name = "Paired")
  cum_comp <- t(apply(ordered_phis, 1, cumsum))
  plot(NA, xlim = xlim, ylim=c(0,1),
       yaxs="i", xaxs="i", 
       # xaxt = "n", yaxt = "n",
       xlab=bquote("Growth rate"~ mu ~ (h^-1)), 
       ylab="Relative proteome composition",
       cex.lab = 1.3)
  for(i in ncol(cum_comp):1){
    col <- colors[i]
    polygon(c(ordered_mu, rev(ordered_mu)), c(rep(0, nrow(cum_comp)), rev(cum_comp[,i])),
            col=col, border=col)
    par(new=TRUE)
  }
  par(new=FALSE)
  # if(legend){
  #   legend('bottomleft', rev(colnames(ordered_phis)), bty = "n",
  #          col = rev(colors), pch = 15, cex=1.1, border = NA)
  # }
  title(dataset)

  # translation sectors
  translation <- get_sector_fraction(c("translation", "protein maturation"),
                                     current_dataset, gene_annot)
  ribosomes <- get_sector_fraction_gene_list(ribo_genes, current_dataset)
  trna_synthetases <- get_sector_fraction_gene_list(trna_syn, current_dataset)
  
  sectors <- list("Translation (w/o tRNA synthetases)" = translation-trna_synthetases,
                  "tRNA synthetases" = trna_synthetases)
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.45), 
               xlim = xlim,
               legend = legend)
  title(dataset)

  # translation sectors
  sectors <- list("Translation (w/o tRNA synthetases)" = translation-trna_synthetases,
                  "Ribosomal proteins" = ribosomes)
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.45), 
               xlim = xlim,
               legend = legend)
  title(dataset)
  
  
  rib_data <- data.frame("dataset"=dataset, "translation"=translation, "ribosomes"=ribosomes)
  if(dataset != "Merged"){
    all_rib_data <- rbind(all_rib_data, rib_data)
  }

  
  # RNA and DNA
  rna <- get_sector_fraction(c("DNA-templated transcription", "RNA processing"),
                             current_dataset, gene_annot)
  dna <- get_sector_fraction(c("DNA biosynthetic process"), current_dataset, gene_annot)
  sectors <- list("Transcription" = rna,
                  "DNA replication" = dna)

  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.1), 
               xlim = xlim,
               legend = legend)
  title(dataset)
  

  # metabolism
  upper <- get_sector_fraction(upper_groups, current_dataset, gene_annot)
  etc <- get_sector_fraction_gene_list(etc_genes, current_dataset)
  #energy <- get_sector_fraction("energy derivation by oxidation of organic compounds", current_dataset, gene_annot)
  sectors <- list("Central Carbon metabolism" = upper,
                  "ETC" = etc)
                  #"TCA+OxOphos" = energy)
                  #"All carbohydrate/energy metabolism" = upper+etc
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.5), 
               xlim = xlim,
               legend = legend)
  title(dataset)


  
  # transport
  gt <- get_sector_fraction(c("D-glucose transmembrane transport"), current_dataset, gene_annot)
  wt <- get_sector_fraction(c("carboxylic acid transport"), current_dataset, gene_annot)

  sectors <- list("Glucose transport" = gt,
                  "Carboxylic acid transport" = wt)
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.06), 
               xlim = xlim,
               legend = legend)
  title(dataset)
  
  
  # AAs, nts, lipids
  lip <- get_sector_fraction(c("lipid biosynthetic process"), current_dataset, gene_annot)
  aas <- get_sector_fraction(c("amino acid biosynthetic process"), current_dataset, gene_annot)
  nts <- get_sector_fraction(c("nucleobase biosynthetic process", "NADPH regeneration"), current_dataset, gene_annot)
  
  sectors <- list("Nucleotide metabolism" = nts,
                  "Amino acid metabolism" = aas,
                  "Lipid metabolism" = lip)
  
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.2), 
               xlim = xlim,
               legend = legend)
  title(dataset)
  
  
  
  # AAs, nts, lipids
  lip <- get_sector_fraction(c("lipid biosynthetic process"), current_dataset, gene_annot)
  sectors <- list("Lipid metabolism" = lip)
  
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.05), 
               xlim = xlim,
               legend = legend)
  title(dataset)
  
  
  
  
  # RNA degradation
  rrnases <- get_sector_fraction(c("rRNA catabolic process"), current_dataset, gene_annot)
  trnases <- get_sector_fraction(c("tRNA decay"), current_dataset, gene_annot)
  mrnases <- get_sector_fraction(c("mRNA catabolic process"), current_dataset, gene_annot)
  
  sectors <- list("rRNases" = rrnases,
                  "tRNases" = trnases,
                  "mRNases"= mrnases)
  plot_and_fit(sectors, 
               current_growth, 
               ylim = c(0,0.008), 
               ylab = "Proteome fraction",
               xlim = xlim,
               bquote("Growth rate"~ mu ~ (h^-1)),
               legend = legend)
  title(dataset)

}

write.csv(all_rib_data, "data/ribosome_data.csv")



dev.off()



library(here)

setwd(here())

ecocyc_annot <- read.table("data/p17_groups_detailed.txt", sep = "\t", header=TRUE)
ecocyc_transporters <- read.table("data/ecoli-transmembrane-transport.txt", sep = "\t", header=TRUE)

p17_groups <- c("translation",
                "lipid biosynthetic process",          
                "amino acid biosynthetic process",
                "nucleobase biosynthetic process",
                "DNA-templated transcription",
                "RNA processing",
                "rRNA catabolic process",                       
                "mRNA catabolic process",                        
                "tRNA decay",
                "DNA biosynthetic process",
                #"generation of precursor metabolites and energy",
                  "glycolytic process",
                  "energy derivation by oxidation of organic compounds",
                  "NADPH regeneration", # PPP
                  "ketone body metabolic process",
                  #"electron transport chain",
                "protein maturation",
                "carboxylic acid transport",
                "D-glucose transmembrane transport")

get_genes <- function(ecocyc_annot, groupname){
  genes <- ecocyc_annot[ecocyc_annot$Gene.Ontology.Terms == groupname, 2]
  genes <- strsplit(genes, " // ")[[1]]
  return(genes)
}
 
# subgroups that will be plotted separately
ribo_genes <- get_genes(ecocyc_annot, "structural constituent of ribosome")
trna_syn <- get_genes(ecocyc_annot, "tRNA aminoacylation")
etc_genes <- get_genes(ecocyc_annot, "electron transport chain") # subgroup of "energy derivation by oxidation of organic compounds",

transporter_genes <- ecocyc_transporters[ecocyc_transporters$New.column == "transmembrane transport", 2]
transporter_genes <- strsplit(transporter_genes, " // ")[[1]]

gene_annot <- list()
for(row in 1:nrow(ecocyc_annot)){
  group <- ecocyc_annot[row,1]
  if(!group %in% p17_groups){
    next
  }
  genes <- strsplit(ecocyc_annot[row,2], " // ")[[1]]
  for(gene in genes){
    
    error_check <- try({gene_annot[[gene]] <- c(gene_annot[[gene]], group)}, silent = TRUE)
    if(class(error_check) == "try-error"){
      gene_annot[[gene]] <- c(group)
    }
  }
}


# add RNase annotation https://pmc.ncbi.nlm.nih.gov/articles/PMC4931109/
for(rrnase_gene in c("b1286", "b4179", "b3164", "b1084")){
  if("rRNA catabolic process" %in% gene_annot[[rrnase_gene]]){next}
  gene_annot[[rrnase_gene]] <- c(gene_annot[[rrnase_gene]], "rRNA catabolic process")
}


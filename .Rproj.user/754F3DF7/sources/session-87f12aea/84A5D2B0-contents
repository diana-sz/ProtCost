library(RColorBrewer)
library(here)
require(dplyr)

setwd(here())

metadata <- read.csv("data/Schmidt_proteomics_table_s23.csv", skip=2, nrows = 26)
metadata <- metadata[metadata$Strain == "BW25113",]
rownames(metadata) <- metadata[,1]

proteomics_raw <- read.csv("data/Schmidt_proteomics_table_s6.csv", skip=1)
start_index <- grep("Mass", colnames(proteomics_raw))
end_index <- grep("Coefficient", colnames(proteomics_raw)) - 1 
colnames(proteomics_raw) <- proteomics_raw[1,]
proteomics_all <- proteomics_raw[2:nrow(proteomics_raw),]

proteomics_masses <- proteomics_all[, start_index:end_index]
proteomics_masses[] <- lapply(proteomics_masses, as.numeric)
proteomics_masses$Bnumber <- proteomics_all$Bnumber

# minimal media
selected_conditions <- c("Glucose", "Acetate", "Fumarate", "Glucosamine", 
                         "Glycerol", "Pyruvate", "Chemostat 0.5", "Chemostat 0.35",
                         "Chemostat 0.20", "Chemostat 0.12", "Xylose", "Mannose",
                         "Galactose ", "Succinate", "Fructose")
schmidt_proteomics <- proteomics_masses[,selected_conditions]

# sum up duplicated rows (same gene)
schmidt_proteomics <- schmidt_proteomics %>% 
  group_by(proteomics_masses$Bnumber) %>% 
  summarise_all(list(sum=sum))

# put into the right format - gene IDs as row names
schmidt_proteomics <- data.frame(schmidt_proteomics)
rownames(schmidt_proteomics) <- schmidt_proteomics$`proteomics_masses.Bnumber`
schmidt_proteomics <- schmidt_proteomics[,2:ncol(schmidt_proteomics)]
colnames(schmidt_proteomics) <- selected_conditions

schmidt_growth_rates <- metadata[colnames(schmidt_proteomics), "Growth.rate..h.1."]


#### Mori data #################################################################
metadata <- read.csv("data/EV3-Samples-2.csv")
c_lim_samples <- metadata[8:22,]
mori_growth_rates <- c_lim_samples$Growth.rate..1.h.

mori_proteomics <- read.csv("data/EV9-AbsoluteMassFractions-2.csv")
mori_proteomics <- mori_proteomics[mori_proteomics$Gene.locus != "", ]
rownames(mori_proteomics) <- mori_proteomics$Gene.locus
mori_proteomics <- mori_proteomics[, c_lim_samples$Sample.ID]



#### Mori data #################################################################
metadata <- read.csv("data/EV2-Samples-1.csv")
target_rows <- c("Acetate", "Carbon 46", "Carbon 61", "Carbon 37",
                 "Carbon 85", "Carbon 50", "Carbon 54", "Carbon 49")
min_media_samples <- metadata[metadata$Short.Description %in% target_rows,]
mori_growth_rates2 <- log(2)/(min_media_samples$Doubling.time..min./60)

mori_proteomics2 <- read.csv("data/EV8-AbsoluteMassFractions-1.csv")
mori_proteomics2 <- mori_proteomics2[mori_proteomics2$Gene.locus != "", ]
rownames(mori_proteomics2) <- mori_proteomics2$Gene.locus
sample_ids <- gsub("-", ".", min_media_samples$Sample.ID)
mori_proteomics2 <- mori_proteomics2[, sample_ids]



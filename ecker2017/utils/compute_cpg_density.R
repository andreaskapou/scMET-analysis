################################################################################
##  Script to computing CpG density per feature in different genomic contexts ##
################################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

#####################
## Define settings ##
#####################
source("../load_settings.R")
io$outdir <- paste0(io$basedir, "stats/features/")

# Genomic contexts
# opts$annos <- "all"
opts$annos <- c(
  "prom_2000_2000",
  "distal_H3K27ac_cortex",
  "H3K4me1_cortex"
)
# Chromosomes
opts$chr <- c("X","Y",1:19)

if (opts$annos == "all") {
  opts$annos <- list.files(io$features, pattern = "\\.bed.gz$") %>%
    gsub(".bed.gz","",.)
}

###########################
## Load feature metadata ##
###########################
anno_dt <- opts$annos %>% map(~ fread(sprintf("%s/%s.bed.gz", io$features,.))) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno")) %>%
  .[chr %in% opts$chr] %>%
  .[, chr := paste0("chr",chr)]

#######################################
## Calculate CpG density per feature ##
#######################################
seq <- getSeq(Mmusculus, anno_dt$chr, anno_dt$start, anno_dt$end + 1)
anno_dt$cpg_density <- dinucleotideFrequency(seq)[,"CG"] / width(seq)

##################
## Save results ##
##################
fwrite(anno_dt, paste0(io$outdir,"/cpg_density_perfeature.txt.gz"),
       col.names = TRUE, quote = FALSE, sep = "\t")

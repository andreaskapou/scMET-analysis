library(data.table)
library(purrr)
library(ggplot2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

#####################
## Define settings ##
#####################

## Define I/O ##
if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/Ecker_2017/settings.R")
  io$outdir <- "/Users/ricard/data/Ecker_2017/mouse/stats/features"
} else if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  source("~/Research/Projects/scMET/code/scMET-analysis/gastrulation/load_settings.R")
  io$outdir <- paste0(io$basedir,"met/stats/features/")
} else {
  stop("Computer not recognised")
}


## Define options ##

# Genomic contexts
# opts$annos <- c(
#   "prom_2000_2000_noncgi",
#   "prom_2000_2000_cgi",
#   "H3K27ac_cortex",
#   "H3K27ac_cortex",
#   "genebody",
#   "LINE",
#   # "LTR",
#   "CGI"
# )
opts$annos <- c("prom_2000_2000", "first_exon")

# Chromosomes
opts$chr <- c("X","Y",1:19)

if (opts$annos == "all") {
  opts$annos <- list.files(io$features, pattern = "\\.bed.gz$") %>% gsub(".bed.gz","",.)
}
# opts$annos <- opts$annos[1:length(opts$anno)-1]


###########################
## Load feature metadata ##
###########################

anno_dt <- opts$annos %>% map(~ fread(sprintf("%s/%s.bed.gz",io$features,.))) %>%
  rbindlist %>% setnames(c("chr","start","end","strand","id","anno")) %>%
  .[chr%in%opts$chr] %>%
  .[,chr:=paste0("chr",chr)]

#######################################
## Calculate CpG density per feature ##
#######################################

# Sanity checks
chr_lengths.dt <- data.table(
  chr = paste0("chr",opts$chr),
  chr_length = seqlengths(Mmusculus) %>% .[paste0("chr",opts$chr)]
)
anno_dt <- merge(anno_dt, chr_lengths.dt, by="chr")

# Filter features that exceed chr  length
anno_dt <- anno_dt[end<chr_length]


# Get sequence
seq <- getSeq(Mmusculus, anno_dt$chr, anno_dt$start, anno_dt$end+1)

# Calculate CpG density
anno_dt$cpg_density <- dinucleotideFrequency(seq)[,"CG"] / width(seq)

##################
## Save results ##
##################

fwrite(anno_dt, paste0(io$outdir,"/cpg_density_perfeature.txt.gz"), col.names=T, quote=F, sep="\t")

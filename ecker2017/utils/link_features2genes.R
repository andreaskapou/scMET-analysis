##########################################
##  Script to linking features to genes ##
##########################################

library(data.table)
library(purrr)
library(ggplot2)

##############
## Settings ##
##############
source("../load_settings.R")
io$outdir <- paste0(io$basedir, "/features/")

# Genomic annotations
opts$annos <- list(
  "prom_2000_2000",
  "distal_H3K27ac_cortex",
  "H3K4me1_cortex",
  "window10000_step10000",
  "window20000_step20000"
)
opts$gene_centric.annotations <- c("prom_2000_2000")
# window length for the overlap
opts$gene_window <- 2e4
opts$gene_window <- 1e5

###############
## Load data ##
###############
# Load gene metadata
gene_metadata <- fread(io$gene.metadata) %>%
  .[,chr := as.factor(sub("chr","",chr))] %>%
  setnames(c("ens_id","symbol"), c("id","gene"))

# Load genomic annotations
feature_metadata <- lapply(opts$annos, function(n)
  fread(sprintf("%s/%s.bed.gz",io$features,n), showProgress = FALSE)
) %>% rbindlist %>% setnames(c("chr","start","end","strand","id","anno"))

################
## Parse data ##
################
# Prepare metadata for the overlap
gene_metadata_filt <- gene_metadata %>%
  .[, c("chr","start","end","gene","id")] %>%
  .[,c("start","end") := list(start - opts$gene_window, end + opts$gene_window)] %>%
  .[, c("chr","start","end","gene")] %>%
  setkey(chr,start,end)

#############
## Overlap ##
#############
# Gene-centric annotations
ov1 <- feature_metadata[anno %in% opts$gene_centric.annotations] %>%
  setnames(c("start","end"), c("feature.start","feature.end")) %>%
  merge(gene_metadata[,c("chr","start","end","id","gene")] %>%
          setnames(c("start","end"), c("gene.start","gene.end")), by = c("chr","id")) %>%
  .[,dist := 0]
# Non gene-centric annotations
ov2 <- foverlaps(
  feature_metadata[!anno %in% opts$gene_centric.annotations] %>% setkey(chr, start, end),
  gene_metadata_filt,
  nomatch = NA
) %>%
  setnames(c("i.start","i.end"),c("feature.start","feature.end")) %>%
  .[,c("gene.start", "gene.end") := list(start + opts$gene_window, end - opts$gene_window)] %>%
  .[,c("start","end") := NULL] %>%
  .[,c("start_dist","end_dist") := list( gene.end - feature.start, gene.start - feature.end)] %>%
  .[,c("start_dist","end_dist") := list( ifelse(end_dist < 0 & start_dist > 0,0,start_dist),
                                         ifelse(end_dist < 0 & start_dist > 0,0,end_dist) )] %>%
  .[,dist := ifelse(abs(start_dist) < abs(end_dist), abs(start_dist), abs(end_dist))] %>%
  .[,c("start_dist","end_dist") := NULL]

# sort columns
cols <- c("chr", "id", "anno", "feature.start", "feature.end", "strand",
          "gene.start", "gene.end", "gene", "dist")
ov1 <- ov1[,cols, with = FALSE]
ov2 <- ov2[,cols, with = FALSE]
ov <- rbind(ov1,ov2)

##########
## Save ##
##########
fwrite(ov, paste0(io$outdir,"/genes2features100K.txt.gz"), quote = FALSE,
       col.names = TRUE, row.names = FALSE, sep = "\t", na = "NA")

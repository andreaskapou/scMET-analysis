#########
## I/O ##
#########
library(parallel)
source("../../load_settings.R")
io$script.boxplot <- "boxplot_hits.R"
io$data_dir <- "~/datasets/scMET_ms/gastrulation/data/"
io$hits_dir <- "~/datasets/scMET_ms/gastrulation/metrna/hits/"
io$out_dir <- "~/datasets/scMET_ms/gastrulation/metrna/analysis/examples/"

#############
## Options ##
#############
opts$anno <- c(
  "prom_2000_2000"
)
# Top number of hits to plot
opts$nhits <- 15
opts$show_hits <- "marker_genes"  # "marker_genes" or "top_hits"

###############
## Load data ##
###############
# Load gene metadata
gene_metadata <- fread(io$gene.metadata) %>% .[,c("ens_id","symbol")]

################
## Parse data ##
################
# Load results from HVF analysis
df <- fread(sprintf("%s/met_hvf.txt.gz", io$hits_dir)) %>%
  as.data.table %>%
  merge(gene_metadata, by = "ens_id")


if (opts$show_hits == "marker_genes") {
  # Marker genes
  marker_genes <- c(
    "Apob", "Cer1", "Cubn", "Dppa2", "Dppa5a", "Fmr1nb", "Id3", "Krt8",
    "Lefty2", "Mesp1", "Morc1", "Spp1", "Trap1a", "Zfp42",
    "Gnas", "Usp29", "Brwd3", "Hs6st2", "Yipf6", "Magt1", "Hcfc1",
    "Pnpt1", "Hnrnph2", "Ar", "Pak1ip1", "Rpl10", "Ubqln2", "Zrsr1")
  df <- df[symbol %in% marker_genes]
  opts$nhits <- 30
} else {
  df <- df %>%
    .[tail_prob > 0.9] %>%
    setorder(-tail_prob, -epsilon)
}
hits <- df$ens_id[1:min(opts$nhits, NROW(df))]

##########
## Plot ##
##########
io$out_dir <- paste0(io$out_dir, "/", opts$show_hits, "/")
if (!dir.exists(io$out_dir)) { dir.create(io$out_dir, recursive = TRUE) }
# Extract all cell lineages
stage_lineage <- unique(sample_metadata$stage_lineage)

mclapply(X = 1:length(hits), function(i) {
#for (i in 1:length(hits)) {
    cmd <- sprintf("Rscript %s --gene %s --met.id %s --met.anno %s --stage_lineage %s --datadir %s --outdir %s",
                   io$script.boxplot, df[i, symbol], df[i, ens_id], opts$anno,
                   paste(stage_lineage, collapse = " "), io$data_dir, io$out_dir)
    system(cmd)
}, mc.cores = 8)

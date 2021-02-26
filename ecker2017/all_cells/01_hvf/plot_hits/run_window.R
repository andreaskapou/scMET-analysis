#########
## I/O ##
#########
source("../../../load_settings.R")
io$script.boxplot <- "boxplot_hits_window.R"
io$marker_genes_file <- "~/datasets/scMET_ms/ecker2017/metadata/ecker2017_marker_genes.csv"
io$hits_dir <- "~/datasets/scMET_ms/ecker2017/all_cells/hvf_window/hits/"
io$out_dir <- "~/datasets/scMET_ms/ecker2017/all_cells/hvf_window/examples/"

#############
## Options ##
#############
opts$anno <- c(
  "window20000_step20000" = "Window 20k"
)
# Top number of hits to plot
opts$nhits <- 15
# Threshold on the tail probability
opts$min_prob <- 0.8
# Choose whether to perform HVF calling using overdispersion \gamma
#  or residual overdispersion \epsilon. The output files will be
#  stored with this prefix.
opts$mode <- "epsilon" # Either "epsilon" or "gamma"
opts$show_hits <- "marker_genes"  # "marker_genes" or "top_hits"

###############
## Load data ##
###############
# Load results from HVF analysis
df <- names(opts$anno) %>%
  map(~ fread(sprintf("%s/hvf_%s_chr1_%s.txt.gz",io$hits_dir, ., opts$mode))
  ) %>% rbindlist

if (opts$show_hits == "marker_genes") {
  df <- df[id %in% c("window20000_step20000_2848",
                     "window20000_step20000_3440")]
  opts$nhits <- 2
}

###################
## Plot HVF hits ##
###################
# Genomic features
dt_sub <- df %>%
  .[tail_prob >= opts$min_prob & is_variable == TRUE] %>%
  setorderv(cols = c("tail_prob", opts$mode), order = -1) %>%
  .[,head(.SD, n = opts$nhits), by = "anno"]


io$out_dir <- paste0(io$out_dir, "/", opts$mode, "/", opts$show_hits, "/")
if (!dir.exists(io$out_dir)) { dir.create(io$out_dir, recursive = TRUE) }

# Boxplots of the DNA methylation rate per cell type
for (i in 1:length(opts$anno)) {
  if (length(dt_sub[anno == names(opts$anno[i]), id]) == 0) { next }
  cmd <- sprintf("Rscript %s --anno %s --id %s --outdir %s", io$script.boxplot,
                 names(opts$anno[i]),
                 paste(dt_sub[anno == names(opts$anno[i]),id],collapse = " "), io$out_dir)
  system(cmd)
}

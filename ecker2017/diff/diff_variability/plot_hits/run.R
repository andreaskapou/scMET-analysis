#########
## I/O ##
#########
source("../../../load_settings.R")
io$script.boxplot <- "boxplot_hits.R"
io$marker_genes_file <- "~/datasets/scMET_ms/ecker2017/metadata/ecker2017_marker_genes.csv"
io$hits_dir <- "~/datasets/scMET_ms/ecker2017/diff/hits/"
io$out_dir <- "~/datasets/scMET_ms/ecker2017/diff/diff_var/examples/"

#############
## Options ##
#############
opts$anno <- c(
  "distal_H3K27ac_cortex" = "Distal H3K27ac"
  #"H3K4me1_cortex" = "H3K4me1",
  #"prom_2000_2000" = "Promoters"
)
# Define groups
opts$groups <- c("Excitatory", "Inhibitory")
# Top number of hits to plot
opts$nhits <- 20
# Threshold on the tail probability
opts$tail_prob_threshold <- 0.8
opts$show_hits <- "top_hits"  # "marker_genes" or "top_hits"
opts$mode <- "epsilon" # Either "epsilon" or "gamma"
# Extract corresponding columns
tail_prob <- list(
  "epsilon" = "epsilon_tail_prob",
  "gamma" = "gamma_tail_prob"
)
test_columns <- list(
  "epsilon" = "epsilon_diff_test",
  "gamma" = "gamma_diff_test"
)
eff_size <- list(
  "epsilon" = "epsilon_change",
  "gamma" = "gamma_LOR"
)
###############
## Load data ##
###############
# Load results from differential analysis
df <- list()
for (an in names(opts$anno)) {
  tmp <- readRDS(sprintf("%s/%s_vb.rds", io$hits_dir, an))
  df[[an]] <- tmp$diff_epsilon_summary %>% as.data.table %>%
    .[, anno_name := opts$anno[[an]]] %>%
    .[, anno := an] %>%
    setnames(c("feature_name"), c("id"))
}
df <- rbindlist(df)

if (opts$show_hits == "marker_genes") {
  # Marker genes
  marker_genes <- fread(io$marker_genes_file)
  # Load mapping between features and genes
  features2genes <- fread(io$features2genes)
  # Merge
  df <- df %>% merge(features2genes[,c("id","anno","gene")], by = c("id","anno")) %>%
    .[gene %in% marker_genes$gene_name] %>%
    .[, .SD[which.max(get(tail_prob[[opts$mode]]))], by = c("gene", "anno")]
  opts$nhits <- 20
}


##########################################################
## Plot hits that highly methylated in Excitatory cells ##
##########################################################
# Genomic features
dt_sub <- df %>%
  .[get(tail_prob[[opts$mode]]) >= opts$tail_prob_threshold &
      get(test_columns[[opts$mode]]) == "Exc+"] %>%
  .[, abs_diff := abs(get(eff_size[[opts$mode]]))] %>%
  setorderv(cols = c(tail_prob[[opts$mode]], "abs_diff"), order = -1) %>%
  .[,head(.SD, n = opts$nhits), by = "anno"]

io$outdir2 <- paste0(io$out_dir, "/", opts$mode, "/", opts$show_hits, "/Exc+/")
if (!dir.exists(io$outdir2)) { dir.create(io$outdir2, recursive = TRUE) }

# Boxplots of the DNA methylation rate per cell type
for (i in 1:length(opts$anno)) {
  if (length(dt_sub[anno == names(opts$anno[i]), id]) == 0) { next }
  cmd <- sprintf("Rscript %s --anno %s --id %s --outdir %s", io$script.boxplot,
                 names(opts$anno[i]),
                 paste(dt_sub[anno == names(opts$anno[i]),id],collapse = " "), io$outdir2)
  system(cmd)
}


##########################################################
## Plot hits that highly methylated in Inhibitory cells ##
##########################################################
# Genomic features
dt_sub <- df %>%
  .[get(tail_prob[[opts$mode]]) >= opts$tail_prob_threshold &
      get(test_columns[[opts$mode]]) == "Inh+"] %>%
  .[, abs_diff := abs(get(eff_size[[opts$mode]]))] %>%
  setorderv(cols = c(tail_prob[[opts$mode]], "abs_diff"), order = -1) %>%
  .[,head(.SD, n = opts$nhits), by = "anno"]

io$outdir2 <- paste0(io$out_dir, "/", opts$mode, "/", opts$show_hits, "/Inh+/")
if (!dir.exists(io$outdir2)) { dir.create(io$outdir2, recursive = TRUE) }

# Boxplots of the DNA methylation rate per cell type
for (i in 1:length(opts$anno)) {
  if (length(dt_sub[anno == names(opts$anno[i]), id]) == 0) { next }
  cmd <- sprintf("Rscript %s --anno %s --id %s --outdir %s", io$script.boxplot,
                 names(opts$anno[i]),
                 paste(dt_sub[anno == names(opts$anno[i]),id],collapse = " "), io$outdir2)
  system(cmd)
}

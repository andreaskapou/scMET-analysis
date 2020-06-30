# Load libraries
suppressPackageStartupMessages(library(scMET))
suppressPackageStartupMessages(library(argparse))

################
## Load settings
source("../../load_settings.R")
# Set to smaller value due to downsampling to small number of cells
opts$min.cells <- 5

################################
## Initialize argument parser ##
p <- ArgumentParser(description = '')
p$add_argument('--anno',       type = "character",      help = 'Genomic context')
p$add_argument('--cells',      type = "integer",        help = 'Number of cells to keep')
p$add_argument('--replicate',  type = "integer",        help = 'Replicate ID (integer)')
p$add_argument('--outdir',     type = "character",      help = 'Output directory')
p$add_argument('--test',       action = "store_true",   help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))
anno <- args$anno
N_cells <- args$cells
replicate <- args$replicate
is_test <- isTRUE(args$test)
outdir <- args$outdir
message("Test mode: ", is_test)

## Define seed for reproducibility between MLE and VB downsampling
set.seed(N_cells + replicate + 10000)
# Subset features (for testing)
if (is_test) { N_feat <- 200 }

####################################
## Load / filter methylation data ##
cat("Loading data...\n")
# Keep only "Inhibitory" cells and then downsample
opts$cells <- sample(sample_metadata[Neuron_type1 == "Inhibitory", ]$sample, N_cells)
Y <- read_filter_ecker_data(filename = sprintf("%s/%s.tsv.gz", io$data_parsed, anno),
                            opts = opts, is_differential = FALSE)
# Subset features in testing mode
if (is_test) {
  print("Test mode activated: subsetting features")
  Y <- Y[Feature %in% head(unique(Y$Feature), n = N_feat)]
}

################################################
## Fit Beta binomial maximum likelihood model ##
print("Fitting beta binomial model using maximum likelihood...")
df <- Y[, bb_mle(cbind(total_reads, met_reads)), by = c("Feature")] %>%
  .[, anno := anno]

##########
## Save ##
rep_dir <- paste0(outdir, "/rep", replicate, "/")
if (!dir.exists(rep_dir)) { dir.create(rep_dir, recursive = TRUE) }
fwrite(df, sprintf("%s/%s_%s_%s_mle.txt.gz", rep_dir, anno, "Inh", N_cells), sep = "\t")
cat("Finished!!\n")

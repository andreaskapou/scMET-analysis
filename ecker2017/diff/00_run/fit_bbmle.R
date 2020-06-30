# Load libraries
suppressPackageStartupMessages(library(scMET))
suppressPackageStartupMessages(library(argparse))

################
## Load settings
source("../../load_settings.R")

################################
## Initialize argument parser ##
################################
p <- ArgumentParser(description = '')
p$add_argument('--anno',    type = "character",    help = 'genomic context')
p$add_argument('--group',   type = "character",    help = 'cell group (Inhibitory or Excitatory)')
p$add_argument('--outdir',  type = "character",    help = 'Output directory')
p$add_argument('--test',    action = "store_true", help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))
anno <- args$anno
group <- args$group
outdir <- args$outdir
is_test <- isTRUE(args$test)
message("Test mode: ", is_test)

################
## Model settings (mostly keeping default parameters)
opts$iter <- ifelse(use_mcmc, 2000, 50000)
opts$algorithm <- ifelse(use_mcmc, "NUTS", "meanfield")
opts$chains <- 1
if (use_mcmc) {
  opts$chains <- ifelse(is_test, 1, 2)
}
opts$n_cores <- opts$chains
# Subset features (for testing)
if (is_test) { N_feat <- 200; opts$iter <- 2000 }

################
# Define broad groups of cells
opts$groups <- unique(sample_metadata$Neuron_type1)


####################################
## Load / filter methylation data ##
cat("Loading data...\n")
Y <- read_filter_ecker_data(filename = sprintf("%s/%s.tsv.gz", io$data_parsed, anno),
                            opts = opts, sample_metadata = sample_metadata,
                            is_differential = TRUE)

####################################
## Keep data from specified group ##
Y <- Y %>% .[Cell %in% sample_metadata[Neuron_type1 == group, sample]]
cat("Keep cells from group ", group, "\n")
cat("Total # cells: ", length(unique(Y$Cell)), "\n")
cat("Total # features: ", length(unique(Y$Feature)), "\n")
# Subset features in testing mode
if (is_test) {
  print("Test mode activated: subsetting features")
  Y <- Y[Feature %in% head(unique(Y$Feature), n = N_feat)]
}

################################################
## Fit Beta binomial maximum likelihood model ##
print("Fitting beta binomial model using maximum likelihood...")
df <- Y[, bb_mle(cbind(total_reads, met_reads)), by = c("Feature")] %>%
  .[, anno := anno] %>% .[, group := group]

##########
## Save ##
##########
fwrite(df, sprintf("%s/%s_%s_mle.txt.gz", outdir, anno, group), sep = "\t")
cat("Finished!!\n")

# Load libraries
suppressPackageStartupMessages(library(scMET))
suppressPackageStartupMessages(library(argparse))
set.seed(123)

################
## Load settings
source("../load_settings.R")

# anno <- "prom_2000_2000"
# outdir <- "~/datasets/scMET_ms/gastrulation/data/"
# is_test <- FALSE
# use_mcmc <- FALSE

################################
## Initialize argument parser ##
p <- ArgumentParser(description = '')
p$add_argument('--anno',       type = "character",      help = 'Genomic context')
p$add_argument('--outdir',     type = "character",      help = 'Output directory')
p$add_argument('--mcmc',       action = "store_true",   help = 'Use MCMC? (default is VB)')
p$add_argument('--test',       action = "store_true",   help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))
anno <- args$anno
is_test <- isTRUE(args$test)
outdir <- args$outdir
use_mcmc <- isTRUE(args$mcmc)
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

####################################
## Load / filter methylation data ##
cat("Loading data...\n")
Y <- read_filter_gastrulation_data(filename = sprintf("%s/%s.tsv.gz", io$met_data_parsed, anno),
                                   opts = opts, is_differential = FALSE)
# Subset features in testing mode
if (is_test) {
  print("Test mode activated: subsetting features")
  Y <- Y[Feature %in% head(unique(met$Feature), n = N_feat)]
}

####################################
# Regression on the mean methylation rate by using the CpG density as covariate.s
X <- read_cpg_density(filename = io$cpg_density, feature_names = unique(Y$Feature),
                      annotation = anno)

######################################
## Fit scMET model
ifelse(use_mcmc, "Fitting model using MCMC...", "Fitting model using VB...")
print(date())
# Fit the model
fit <- scmet(Y = Y, X = X, L = 4, use_mcmc = use_mcmc, use_eb = TRUE,
             iter = opts$iter, algorithm = opts$algorithm, output_samples = 2000,
             chains = opts$chains, s_wmu = 2, s_mu = 1.5, s_wgamma = 2, a_sgamma = 2,
             b_sgamma = 3, rbf_c = 1, init_using_eb = TRUE, tol_rel_obj = 1e-04,
             n_cores = opts$n_cores, lambda = 4)
print(date())

##########
## Save ##
cat("Storing results ...\n")
# Extract parameter summaries and store them in separate txt file
df <- extract_param_summaries(scmet_obj = fit, anno = anno)
mode <- ifelse(use_mcmc, "mcmc", "vb")
if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
saveRDS(fit, file = sprintf("%s/%s_%s.rds", outdir, anno, mode))
fwrite(df, sprintf("%s/%s_%s.txt.gz", outdir, anno, mode), sep = "\t")
cat("Finished!!\n")

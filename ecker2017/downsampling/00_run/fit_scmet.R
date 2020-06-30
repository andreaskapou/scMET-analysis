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
p$add_argument('--mcmc',       action = "store_true",   help = 'Use MCMC? (default is VB)')
p$add_argument('--test',       action = "store_true",   help = 'Testing mode')
args <- p$parse_args(commandArgs(TRUE))
anno <- args$anno
N_cells <- args$cells
replicate <- args$replicate
is_test <- isTRUE(args$test)
outdir <- args$outdir
use_mcmc <- isTRUE(args$mcmc)
message("Test mode: ", is_test)

## Define seed for reproducibility between MLE and VB downsampling
set.seed(N_cells + replicate + 10000)

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
# Keep only "Inhibitory" cells and then downsample
opts$cells <- sample(sample_metadata[Neuron_type1 == "Inhibitory", ]$sample, N_cells)
Y <- read_filter_ecker_data(filename = sprintf("%s/%s.tsv.gz", io$data_parsed, anno),
                            opts = opts, is_differential = FALSE)
# Subset features in testing mode
if (is_test) {
  print("Test mode activated: subsetting features")
  Y <- Y[Feature %in% head(unique(Y$Feature), n = N_feat)]
}

####################################
# Regression on the mean methylation rate by using the CpG density as covariate.s
X <- read_cpg_density(filename = io$cpg_density, feature_names = unique(Y$Feature),
                      annotation = anno)

######################################
## Fit scMET model ##
ifelse(use_mcmc, "Fitting model using MCMC...", "Fitting model using VB...")
print(date())
# Fit the model
fit <- scmet(Y = Y, X = X, L = 4, use_mcmc = use_mcmc, use_eb = TRUE,
             iter = opts$iter, algorithm = opts$algorithm, output_samples = 2000,
             chains = opts$chains, s_wmu = 2, s_mu = 1.5, s_wgamma = 2, a_sgamma = 2,
             b_sgamma = 3, rbf_c = 1, init_using_eb = TRUE, tol_rel_obj = 1e-04,
             n_cores = opts$n_cores, lambda = 4)
print(date())


################################################################
## Extract summary statistics for the posterior distributions ##
# Function for computing posterior summary
posterior_summary <- function(dt){
  return(data.table("posterior_median" = matrixStats::colMedians(dt),
                    "posterior_sd" = apply(dt, 2, sd)))
}
dt_mu <- posterior_summary(fit$posterior$mu) %>% setnames(c("mu_median","mu_sd"))
dt_gamma <- posterior_summary(fit$posterior$gamma) %>% setnames(c("gamma_median","gamma_sd"))
dt_epsilon <- posterior_summary(fit$posterior$epsilon) %>% setnames(c("epsilon_median","epsilon_sd"))
summ_stats <- fit$Y[, list(gauss_mean = mean(met_reads / total_reads),
                           gauss_var = var(met_reads / total_reads),
                           cpgs = mean(total_reads), cells = .N), by = c("Feature")]
df <- do.call("cbind", list(dt_mu, dt_gamma, dt_epsilon, summ_stats)) %>% as.data.table %>%
  .[, Feature := fit$feature_names] %>% .[, anno := anno]


##########
## Save ##
cat("Storing results ...\n")
mode <- ifelse(use_mcmc, "mcmc", "vb")
rep_dir <- paste0(outdir, "/rep", replicate, "/")
if (!dir.exists(rep_dir)) { dir.create(rep_dir, recursive = TRUE) }
saveRDS(fit, file = sprintf("%s/%s_%s_%s_%s.rds", rep_dir, anno, "Inh", N_cells, mode))
fwrite(df, sprintf("%s/%s_%s_%s_%s.txt.gz", rep_dir, anno, "Inh", N_cells, mode), sep = "\t")
cat("Finished!!\n")

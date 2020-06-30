##------------------------
# Fit scMET model for synthetic data.
##------------------------
suppressPackageStartupMessages(library(scMET))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(microbenchmark))

##------------------
# I/O
##------------------
io <- list()
if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  io$out_dir <- "~/datasets/scMET_ms/synthetic/scalability/data/"
} else if (grepl("ecdf.ed.ac.uk", Sys.info()['nodename'])) {
  io$out_dir <- "/exports/igmm/eddie/ckapoura-XDF/scMET_ms/synthetic/scalability/data/"
} else {
  stop("Computer not recognised")
}

##------------------------
# Initialise argument parser
##------------------------
p <- ArgumentParser(description = '')
p$add_argument('--replicate',   type = "integer",         help = 'Replicate ID (integer)')
p$add_argument('--cells',       type = "integer",         help = 'Number of cells (integer)')
p$add_argument('--features',    type = "integer",         help = 'Number of features (integer)')
p$add_argument('--mcmc',        action = "store_true",    help = 'Use MCMC? (default is VB)')
p$add_argument('--test',        action = "store_true",    help = 'Testing mode?')
#p$add_argument('--seed',       type = "integer",         help = 'Random seed')

# Parse arguments
args <- p$parse_args(commandArgs(TRUE))
r <- args$replicate     # Number of replicates
N_cells <- args$cells   # Number of cells
N_feat <- args$features # Number of features
use_mcmc <- isTRUE(args$mcmc)

message("Use MCMC: ", use_mcmc)
message("Replicate: ", r)
message("Number of cells: ", N_cells)
message("Number of features: ", N_feat)


##------------------------
# Initialize prior parameters (mostly set to default values)
##------------------------
L <- 4
use_eb <- TRUE
iter <- ifelse(use_mcmc, 3000, 20000)
algorithm <- ifelse(use_mcmc, "NUTS", "meanfield")
output_samples <- 1000
chains <- ifelse(use_mcmc, 1, 1)
s_wmu <- 2
s_mu <- 1.5
s_wgamma <- 2
a_sgamma <- 2
b_sgamma <- 3
rbf_c <- 1
init_using_eb <- TRUE
tol_rel_obj <- 1e-4
n_cores <- chains
lambda <- 4

# Testing mode
if (isTRUE(args$test)) {
  message("Testing mode...")
  r <- 1
  N_cells <- 20
  N_feat <- 100
  iter <- 2000
}


##------------------------
# Load simulated data
##------------------------
cat("Loading simulated data...\n")
N_cpgs <- 15  # Number of CpGs
# Set replication directory
rep_dir <- paste0(io$out_dir, "rep", r, "/")
sim_dt <- readRDS(paste0(rep_dir, "data_feat", N_feat,
                         "_cells", N_cells, "_cpgs", N_cpgs, ".rds"))


##------------------------
# Fit model to data
##------------------------
cat("Performing inference...\n")
# Start time
print(date())
seed <- sample.int(.Machine$integer.max, 1)
ptm <- proc.time()
microbenchmark_time_est <- microbenchmark(
  scmet(Y = sim_dt$Y, X = sim_dt$X, L = L, use_mcmc = use_mcmc, use_eb = use_eb,
        iter = iter, algorithm = algorithm, output_samples = output_samples,
        chains = chains, s_wmu = s_wmu, s_mu = s_mu,
        s_wgamma = s_wgamma, a_sgamma = a_sgamma, b_sgamma = b_sgamma,
        rbf_c = rbf_c, init_using_eb = init_using_eb,
        tol_rel_obj = tol_rel_obj, n_cores = n_cores,
        lambda = lambda, seed = seed), times = 1L)
proc_time_est <- proc.time() - ptm
print(date())

##------------------------
# Store results locally
##------------------------
cat("Storing results locally")
obj <- list(microbenchmark_time_est = microbenchmark_time_est,
            proc_time_est = proc_time_est)
saveRDS(object = obj, file = paste0(rep_dir, "scmet_feat", N_feat,
                                    "_cells", N_cells, "_cpgs", N_cpgs,
                                    "_mcmc", use_mcmc, "_iter", iter, ".rds"))

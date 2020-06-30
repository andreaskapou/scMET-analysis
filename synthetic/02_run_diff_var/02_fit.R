##------------------------
# Fit scMET model for synthetic data.
##------------------------
suppressPackageStartupMessages(library(scMET))
suppressPackageStartupMessages(library(argparse))

##------------------
# I/O
##------------------
io <- list()
if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  io$out_dir <- "~/datasets/scMET/synthetic/diff_var/"
} else if (grepl("ecdf.ed.ac.uk", Sys.info()['nodename'])) {
  io$out_dir <- "/exports/igmm/eddie/ckapoura-XDF/scMET-analysis/data/synthetic/diff_var/"
} else {
  stop("Computer not recognised")
}

##------------------------
# Initialise argument parser
##------------------------
p <- ArgumentParser(description = '')
p$add_argument('--replicate',   type = "integer",         help = 'Replicate ID (integer)')
p$add_argument('--cells',       type = "integer",         help = 'Number of cells (integer)')
p$add_argument('--cpgs',        type = "integer",         help = 'Number of CpGs (integer)')
p$add_argument('--oddsratio',   type = "double",          help = 'Odds Ratio threshold (double)')
p$add_argument('--mcmc',        action = "store_true",    help = 'Use MCMC? (default is VB)')
p$add_argument('--test',        action = "store_true",    help = 'Testing mode?')

# Parse arguments
args <- p$parse_args(commandArgs(TRUE))
r <- args$replicate    # Number of replicates
N_cells <- args$cells  # Number of cells
N_cpgs <- args$cpgs    # Number of CpGs
OR_change_gamma <- args$oddsratio # Odds Ratio threshold
use_mcmc <- isTRUE(args$mcmc)

message("Use MCMC: ", use_mcmc)
message("Replicate: ", r)
message("Number of cells: ", N_cells)
message("Number of CpGs: ", N_cpgs)
message("Odds Ratio threshold: ", OR_change_gamma)

##------------------------
# Initialize prior parameters (mostly set to default values)
##------------------------
L <- 4
use_eb <- TRUE
iter <- ifelse(use_mcmc, 2000, 20000)
algorithm <- ifelse(use_mcmc, "NUTS", "meanfield")
output_samples <- 2000
chains <- ifelse(use_mcmc, 2, 1)
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

#args$test <- TRUE
# Testing mode
if (isTRUE(args$test)) {
  message("Testing mode...")
  r <- 1
  OR_change_gamma <- 5
  N_cells <- 20
  iter <- 2000
}

##------------------------
# Load simulated data
##------------------------
cat("Loading simulated data...\n")
N_feat <- 300 # Number of features
rep_dir <- paste0(io$out_dir, "rep", r, "/")
sim_dt <- readRDS(paste0(rep_dir, "data_ORgamma", OR_change_gamma, "_feat", N_feat,
                         "_cells", N_cells, "_cpgs", N_cpgs, ".rds"))

##------------------------
# Fit group A
##------------------------
cat("Inference for group A\n")
print(date())
seed_A <- sample.int(.Machine$integer.max, 1)
scmet_A <- scmet(Y = sim_dt$scmet_dt_A$Y, X = sim_dt$scmet_dt_A$X,
                 L = L, use_mcmc = use_mcmc, use_eb = use_eb, iter = iter,
                 algorithm = algorithm, output_samples = output_samples,
                 chains = chains, s_wmu = s_wmu, s_mu = s_mu,
                 s_wgamma = s_wgamma, a_sgamma = a_sgamma, b_sgamma = b_sgamma,
                 rbf_c = rbf_c, init_using_eb = init_using_eb,
                 tol_rel_obj = tol_rel_obj, n_cores = n_cores,
                 lambda = lambda, seed = seed_A)
print(date())

##------------------------
# Fit group B
##------------------------
cat("Inference for group B")
print(date())
seed_B <- sample.int(.Machine$integer.max, 1)
scmet_B <- scmet(Y = sim_dt$scmet_dt_B$Y, X = sim_dt$scmet_dt_B$X,
                 L = L, use_mcmc = use_mcmc, use_eb = use_eb, iter = iter,
                 algorithm = algorithm, output_samples = output_samples,
                 chains = chains, s_wmu = s_wmu, s_mu = s_mu,
                 s_wgamma = s_wgamma, a_sgamma = a_sgamma, b_sgamma = b_sgamma,
                 rbf_c = rbf_c, init_using_eb = init_using_eb,
                 tol_rel_obj = tol_rel_obj, n_cores = n_cores,
                 lambda = lambda, seed = seed_B)
print(date())

##------------------------
# Store results locally
##------------------------
# Set model mode
model <- ifelse(L == 1, "_Noreg", "")
cat("Storing results locally")
obj <- list(scmet_A = scmet_A, scmet_B = scmet_B, sim_dt = sim_dt)
saveRDS(object = obj, file = paste0(rep_dir, "/scmet_ORgamma", OR_change_gamma,
                                    "_feat", N_feat, "_cells", N_cells,
                                    "_cpgs", N_cpgs, "_mcmc", use_mcmc, model, ".rds"))

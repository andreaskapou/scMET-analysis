##------------------------
# Simulate synthetic data
# NOTE: There is no downstream analysis performed in this script.
##------------------------
suppressPackageStartupMessages(library(scMET))

##------------------
# I/O
##------------------
io <- list()
if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  io$out_dir <- "~/datasets/scMET/synthetic/diff_var/"
} else if (grepl("ecdf.ed.ac.uk", Sys.info()['nodename'])) {
  io$out_dir <- "/exports/igmm/eddie/ckapoura-XDF/scMET_analysis/data/synthetic/diff_var/"
} else {
  stop("Computer not recognised")
}

##------------------------
# Simulate methylation data from both groups
##------------------------
w_mu <- c(-0.5, -1.5)
s_mu <- 1
w_gamma <- c(-1.2, -.3, 1.1, -.9)
s_gamma <- 0.25
L <- 4
diff_feat_prcg_gamma <- 0.15
OR_change_gamma <- c(2, 3, 5)
N_feat <- 300
N_cpgs <- c(15, 50)
N_cells <- c(20, 50, 100, 200, 500, 1000)
total_rep <- 10

for (r in 1:total_rep) {
  # Set replication directory
  rep_dir <- paste0(io$out_dir, "rep", r, "/")
  cat("Directory ", rep_dir, "\n")
  if (!dir.exists(rep_dir)) { dir.create(rep_dir, recursive = TRUE) }
  for (or in OR_change_gamma) {
    for (c in N_cells) {
      for (cpg in N_cpgs) {
        for (f in N_feat) {
          cat("Running simulation\n")
          seed <- sample.int(.Machine$integer.max, 1)
          sim_dt <- scmet_simulate_diff(N_feat = f, N_cells = c, N_cpgs = cpg, L = L,
                                        diff_feat_prcg_mu = 0,
                                        diff_feat_prcg_gamma = diff_feat_prcg_gamma,
                                        OR_change_gamma = or,
                                        X = NULL, w_mu = w_mu, s_mu = s_mu,
                                        w_gamma = w_gamma, s_gamma = s_gamma,
                                        rbf_c = 1, cells_range = c(0.4, 0.8),
                                        cpgs_range = c(0.4, 0.8), seed = seed)
          ##------------------------
          # Store results locally
          ##------------------------
          saveRDS(object = sim_dt, file = paste0(rep_dir, "/data_ORgamma", or, "_feat", f,
                                                 "_cells", c, "_cpgs", cpg, ".rds"))
        }
      }
    }
  }
}

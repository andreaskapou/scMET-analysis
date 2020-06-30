##------------------------
# Simulate synthetic data
##------------------------
suppressPackageStartupMessages(library(scMET))

##------------------
# I/O
##------------------
io <- list()
if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  io$out_dir <- "~/datasets/scMET/synthetic/scalability/data/"
} else if (grepl("ecdf.ed.ac.uk", Sys.info()['nodename'])) {
  io$out_dir <- "/exports/igmm/eddie/ckapoura-XDF/scMET_ms/synthetic/scalability/data/"
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
N_cpgs <- 15



#N_cells <- c(20, 50, 100, 200, 500, 1000, 2000, 5000)
N_cells <- c(200)
#N_feat <- c(500)
N_feat <- c(20, 50, 100, 200, 1000, 2000, 5000)

total_rep <- 5
for (r in 1:total_rep) {
  # Set replication directory
  rep_dir <- paste0(io$out_dir, "rep", r, "/")
  cat("Directory ", rep_dir, "\n")
  if (!dir.exists(rep_dir)) { dir.create(rep_dir, recursive = TRUE) }
  for (c in N_cells) {
    for (f in N_feat) {
      cat("Running simulation\n")
      seed <- sample.int(.Machine$integer.max, 1)
      sim_dt <- scmet_simulate(N_feat = f, N_cells = c, N_cpgs = N_cpgs,
                               L = L, X = NULL, w_mu = w_mu, s_mu = s_mu,
                               w_gamma = w_gamma, s_gamma = s_gamma,
                               rbf_c = 1, cells_range = c(0.4, 0.8),
                               cpgs_range = c(0.4, 0.8), seed = seed)
      ##------------------------
      # Store results locally
      ##------------------------
      saveRDS(object = sim_dt, file = paste0(rep_dir, "/data_feat", f,
                                             "_cells", c, "_cpgs", N_cpgs, ".rds"))
    }
  }
}

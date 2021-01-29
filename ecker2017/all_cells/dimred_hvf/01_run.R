# Load libraries
library(parallel)

# Source settings
source("../../load_settings.R")
out_dir <- "/Users/ckapoura/datasets/scMET_ms/ecker2017/all_cells/dimred_hvf/data/"
if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }

#############
## Options ##
# Define genomic contexts
opts$anno <- c(
  "prom_2000_2000",
  "distal_H3K27ac_cortex",
  "H3K4me1_cortex"
)

# Define which model to fit
opts$models <- c(
  "scmet",
  "binomial",
  "gaussian",
  "normdisp",
  "random"
)

# Number of highly variable features
opts$number.hvf <- c(seq(50, 1000, by = 50), 2000)

#########
## Run ##
mclapply(X = opts$models, function(model) {
  mclapply(X = opts$anno, function(an) {
    mclapply(X = opts$number.hvf, FUN = function(hvf) {
      outprefix <- sprintf("%s/%s_%s_%d", out_dir, model, an, hvf)
      cmd <- sprintf("Rscript 02_dimred_cluster.R --model %s --anno %s --hvf %d --outprefix %s",
                     model, an, hvf, outprefix)
      system(cmd)
    }, mc.cores = 3)
  }, mc.cores = 3)
}, mc.cores = 3)

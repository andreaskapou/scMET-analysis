suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(BASiCS))
suppressPackageStartupMessages(library(SingleCellExperiment))
set.seed(123)

#####################
## Define settings ##
#####################
source("../load_settings.R")
outdir <- "/Users/ckapoura/datasets/scMET_ms/gastrulation/data/"
if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
table(sample_metadata$stage)

###############
## Load data ##
###############
# SingleCellExperiment object
sce <- readRDS(io$rna)[, sample_metadata$id_rna]

# Keep high quality genes
gex <- BASiCS_Filter(counts(sce),
                     MinTotalCountsPerCell = 10,
                     MinTotalCountsPerGene = 10,
                     MinCellsWithExpression = 10,
                     MinAvCountsPerCellsWithExpression = 5)
# Add a single batch
batch <- rep(1, NCOL(sce))
gex <- newBASiCS_Data(gex$Counts, BatchInfo = batch)

# Run BASiCS to model scRNA-seq data
basics_obj <- BASiCS_MCMC(
  Data = gex,
  N = 20000,
  Thin = 10,
  Burn = 10000,
  Regression = TRUE,
  PriorParam = BASiCS_PriorParam(gex, PriorMu = "EmpiricalBayes"),
  WithSpikes = FALSE
  #Threads = options()[["mc.cores"]]
)
# Save BASiCS object
saveRDS(basics_obj, file = paste0(outdir,"/rna_basics.rds"))

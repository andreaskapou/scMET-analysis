# Load libraries
suppressPackageStartupMessages(library(scMET))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(matrixStats))


# Read CpG density information
# will be used as covariate for regression on mean methylation.
read_cpg_density <- function(filename, feature_names, annotation) {
  cpg_dens <- fread(filename) %>% .[anno == annotation] %>%
    .[, c("id", "anno", "cpg_density")] %>% setnames("id", "Feature") %>%
    .[Feature %in% feature_names] %>%
    .[, logit_cpg := qlogis(cpg_density + 1e-3)] %>%
    .[, logit_cpg_centred := scale(logit_cpg, scale = FALSE), by = c("anno")] %>%
    setorder(Feature)
  # Create covariates X, with intercept 1 and logit CpG density
  X <- cbind(rep(1, NROW(cpg_dens)), cpg_dens$logit_cpg_centred)
  colnames(X) <- c("intercept", "cpg_density_centred")
  rownames(X) <- cpg_dens$Feature
  return(X)
}

# Read methylation data and perform filtering
read_filter_gastrulation_data <- function(filename, opts, sample_metadata = NULL,
                                          is_differential = FALSE) {
  if (length(opts$groups) > 2) {
    stop("Currently this performs analysis for at most 2 groups!")
  }
  if (is.null(sample_metadata) & is_differential) {
    stop("Need to provide sample metadata object for differential analysis!")
  }
  # What grouping we will perform
  if (is_differential) {
    met_groups <- c("Feature", "anno", "group")
  } else {
    met_groups <- c("Feature", "anno")
  }

  #############################
  ## Read methylation data   ##
  #############################
  met <- fread(filename) %>% .[V1 %in% opts$cells] %>%
    setnames(c("Cell", "Feature", "anno", "met_reads", "total_reads", "rate"))
  cat("Total # cells: ", length(unique(met$Cell)), "\n")
  cat("Total # features: ", length(unique(met$Feature)), "\n")

  # If differential testing, add group information to met object
  if (is_differential) {
    # Define broad groups of cells
    met <- met %>% .[, group := opts$groups[1]]
    for (m in 2:length(opts$groups)) {
      cells <- sample_metadata[stage == opts$groups[m], id_met]
      met <- met %>% .[Cell %in% cells, group := opts$groups[m] ]
    }
  }

  #############################
  ## Filter methylation data ##
  #############################
  cat("Filtering data...\n")
  cat("By number of CpGs and minimum number of cells...\n")
  # Filter features by number of CpGs
  met <- met[total_reads >= opts$min.cpgs]
  # Filter features by minimum number of cells
  met <- met[,Ncells := .N, by = met_groups] %>%
    .[Ncells >= opts$min.cells] %>% .[,Ncells := NULL]
  cat("Total # cells: ", length(unique(met$Cell)), "\n")
  cat("Total # features: ", length(unique(met$Feature)), "\n")

  # Keep regions that are not in methylation extremes
  cat("Removing highly/lowly mean methylated features...\n")
  met <- met[, mean_rate := mean(rate/100), by = c("Feature", "anno")] %>%
    .[mean_rate > opts$met_rate_low & mean_rate < opts$met_rate_high] %>%
    .[, mean_rate := NULL]
  cat("Total # cells: ", length(unique(met$Cell)), "\n")
  cat("Total # features: ", length(unique(met$Feature)), "\n")

  # Remove non-variable features
  cat("Removing non-variable features below", opts$var.thresh,  "...\n")
  met <- met[, var := var(rate/100), by = c("Feature", "anno")] %>%
    .[, .SD[var > opts$var.thresh], by = "anno"] %>% .[,var := NULL]
  cat("Total # cells: ", length(unique(met$Cell)), "\n")
  cat("Total # features: ", length(unique(met$Feature)), "\n")

  ## Remove again features by minimum number of cells
  # Filter features by minimum number of cells
  cat("By minimum number of cells (again) ...\n")
  met <- met[,Ncells := .N, by = met_groups] %>%
    .[Ncells >= opts$min.cells] %>% .[,Ncells := NULL]
  cat("Total # cells: ", length(unique(met$Cell)), "\n")
  cat("Total # features: ", length(unique(met$Feature)), "\n")

  if (is_differential) {
    # Obtain intersection of features across groups
    cat("Keeping features that intersect across groups ...\n")
    feature_ids <- fintersect(met[group == opts$groups[1], "Feature", with = FALSE],
                              met[group == opts$groups[2], "Feature", with = FALSE] )
    met <- met[Feature %in% feature_ids$Feature, ]
    cat("Total # cells: ", length(unique(met$Cell)), "\n")
    cat("Total # features: ", length(unique(met$Feature)), "\n")
  }

  ############################
  ## Parse methylation data ##
  ############################
  met <- met[, c("Feature", "Cell", "total_reads", "met_reads")]
  # Reorder data by feature names.
  met <- met[order(Feature), ]

  return(met)
}

# Extract parameter summaries
extract_param_summaries <- function(scmet_obj, anno) {
  # Function for computing posterior summary
  posterior_summary <- function(dt){
    return(data.table("posterior_median" = colMedians(dt),
                      "posterior_sd" = apply(dt, 2, sd)))
  }
  dt_mu <- posterior_summary(scmet_obj$posterior$mu) %>%
    setnames(c("mu_median","mu_sd"))
  dt_gamma <- posterior_summary(scmet_obj$posterior$gamma) %>%
    setnames(c("gamma_median","gamma_sd"))
  dt_epsilon <- posterior_summary(scmet_obj$posterior$epsilon) %>%
    setnames(c("epsilon_median","epsilon_sd"))
  summ_stats <- scmet_obj$Y[, list(gauss_mean = mean(met_reads / total_reads),
                                   gauss_var = var(met_reads / total_reads),
                                   binom_mean = mean(met_reads / total_reads),
                                   binom_var = mean(met_reads / total_reads) *
                                     (1 - mean(met_reads / total_reads)),
                                   cpgs = mean(total_reads), cells = .N),
                            by = c("Feature")]
  df <- do.call("cbind", list(dt_mu, dt_gamma, dt_epsilon, summ_stats)) %>%
    as.data.table %>% .[, Feature := scmet_obj$feature_names] %>% .[, anno := anno]
  return(df)
}

################
## Define I/O ##
################
io <- list()
if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ckapoura/datasets/gastrulation/"
  io$gene.metadata <- "~/datasets/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else if (grepl("ecdf.ed.ac.uk", Sys.info()['nodename'])) {
  io$basedir <- "/exports/igmm/eddie/ckapoura-XDF/gastrulation/"
} else{
  stop("Computer not recognised")
}
io$metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met_data_parsed <- paste0(io$basedir,"/met/parsed")
io$rna <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")
io$features  <- paste0(io$basedir, "/features/filt")
io$cpg_density  <- paste0(io$basedir, "/met/stats/features/cpg_density_perfeature.txt.gz")

####################
## Define options ##
####################
opts <- list()

# Data filtering options
opts$min.cells <- 15          # Minimum # of cells with CpG coverage
opts$min.cpgs <- 3            # Minimum # of CpGs per feature
opts$var.thresh <- 1e-4       # Keep features that have variance above % threshold
opts$met_rate_low <- 0.1      # Threshold for minimum mean methylation rate
opts$met_rate_high <- 0.9     # Threshold for maximum mean methylation rate

# Filtering lineages
opts$remove_lineage_10x_2 <- c(
  "Visceral_endoderm"
)

# ##########################
# ## Load sample metadata ##
# ##########################
sample_metadata <- fread(io$metadata) %>%
  .[,stage_lineage := paste(stage, lineage10x_2,sep = "_")] %>%
  .[!(lineage10x_2 %in% opts$remove_lineage_10x_2)] %>%
  .[pass_metQC == TRUE & pass_rnaQC == TRUE] %>%
  .[,c("sample", "id_met", "id_rna", "embryo", "plate", "stage",
       "lineage10x_2", "stage_lineage")] %>%
  na.omit()

# Define which cells to use
opts$cells <- sample_metadata$id_met

# fwrite(sample_metadata,
#        file = "~/datasets/scMET_ms/gastrulation/metadata/Table_S4_Gastrulation_sample_metadata.csv",
#        sep = "\t")

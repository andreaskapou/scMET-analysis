# Load libraries
suppressPackageStartupMessages(library(scMET))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(magrittr))

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
read_filter_ecker_data <- function(filename, opts, sample_metadata = NULL,
                                   is_differential = FALSE) {
  if (length(opts$groups) > 2) {
    stop("Currently we perform analysis for 2 groups!")
  }
  if (is.null(sample_metadata) & is_differential) {
    stop("Need to provide sample metadata object for differential analysis!")
  }
  if (is_differential & is_downsampling) {
    stop("Cannot perform both differential and downsampling analysis!")
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
    setnames(c("Cell", "Feature", "anno", "rate", "total_reads")) %>%
    .[, met_reads := round((rate/100) * total_reads)]
  cat("Total # cells: ", length(unique(met$Cell)), "\n")
  cat("Total # features: ", length(unique(met$Feature)), "\n")

  # If differential testing, add group information to met object
  if (is_differential) {
    met <- met[, group := opts$groups[1]]
    for (m in 2:length(opts$groups)) {
      cells <- sample_metadata[Neuron_type1 == opts$groups[m], sample]
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

  # Keep regions that are not in mean methylation extremes
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


################
## Define I/O ##
################
io <- list()
if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ckapoura/datasets/ecker2017/mouse/"
} else if (grepl("ecdf.ed.ac.uk", Sys.info()['nodename'])) {
  io$basedir <- "/exports/igmm/eddie/ckapoura-XDF/ecker2017/mouse"
} else{
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$data_parsed <- paste0(io$basedir,"/parsed_CG")
io$features  <- paste0(io$basedir, "/features/filt")
io$cpg_density  <- paste0(io$basedir, "/stats/features/cpg_density_perfeature.txt.gz")

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

# Define excitatory and inhibitory cell types
opts$excitatory <- c("mL2/3","mL4","mL5-1","mDL-1","mDL-2","mDL-3","mL6-1","mL6-2","mL5-2")
opts$inhibitory <- c("mNdnf-1","mNdnf-2","mVip","mPv","mSst-1","mSst-2")
opts$neurons <- c(opts$excitatory, opts$inhibitory)

# Define which cells to use
opts$cells <- fread(io$metadata) %>% .[specie == "Mus_musculus", sample]

###################
## Define colors ##
###################

# Colors for level 1 cell types
opts$colors1 <- c(
  "Excitatory" = "#e31a1c",
  "Inhibitory" = "#1f78b4"
)

# Colors for level 3 cell types
opts$colors3 <- c(
  "mL5-1" = "#696969",  "mL5-2" = "#FEB24C",
  "mL6-1" = "#FED976",  "mL6-2" = "#FFEDA0",
  "mDL-12" = "#FC4E2A", "mDL-3" = "#FD8D3C",
  "mL2/3" = "#E31A1C",  "mL4" = "#8B4513",
  "mVip" = "#8c96c6",   "mSst-12" = "#bae4bc",
  "mPv" = "#7bccc4",    "mNdnf-12" = "#2b8cbe"
)

##########################
## Load sample metadata ##
##########################
sample_metadata <- fread(io$metadata) %>%
  .[`Neuron type` %in% opts$neurons] %>%
  .[,c("sample","specie","Neuron type","Laminar layer", "mCG/CG", "mCH/CH", "Coverage (%)")] %>%
  .[specie == "Mus_musculus"] %>%
  setnames(c("Neuron type", "Laminar layer", "Coverage (%)"), c("Neuron_type", "Laminar_layer", "Coverage(%)"))

# Level 1 cell types
sample_metadata[,"Neuron_type1" := ifelse(Neuron_type %in% opts$excitatory, "Excitatory", "Inhibitory")]

# Level 2 cell types
sample_metadata %>% .[,"Neuron_type2" := as.character(Neuron_type)] %>%
  .[,"Neuron_type2" := ifelse(Neuron_type2 %in% c("mDL-1","mDL-2","mDL-3"), "mDL-123", Neuron_type2)] %>%
  .[,"Neuron_type2" := ifelse(Neuron_type2 %in% c("mNdnf-1","mNdnf-2"), "mNdnf-12", Neuron_type2)] %>%
  .[,"Neuron_type2" := ifelse(Neuron_type2 %in% c("mSst-1","mSst-2"), "mSst-12", Neuron_type2)] %>%
  .[,"Neuron_type2" := ifelse(Neuron_type2 %in% c("mL6-1","mL6-2"), "mL6-12", Neuron_type2)] %>%
  .[,"Neuron_type2" := ifelse(Neuron_type2 %in% c("mL5-1","mL5-2"), "mL5-12", Neuron_type2)]

# Level 3 cell types
sample_metadata %>% .[,"Neuron_type3" := as.character(Neuron_type)] %>%
  .[,"Neuron_type3" := ifelse(Neuron_type3 %in% c("mDL-1", "mDL-2"), "mDL-12", Neuron_type3)]  %>%
  .[,"Neuron_type3" := ifelse(Neuron_type3 %in% c("mNdnf-1", "mNdnf-2"), "mNdnf-12", Neuron_type3)] %>%
  .[,"Neuron_type3" := ifelse(Neuron_type3 %in% c("mSst-1", "mSst-2"), "mSst-12", Neuron_type3)]


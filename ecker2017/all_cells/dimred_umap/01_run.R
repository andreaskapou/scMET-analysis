# Load libraries
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(rhdf5))
set.seed(12345)

# Source settings and I/O
source("../../load_settings.R")
io$fitdir <- "~/datasets/scMET_ms/ecker2017/all_cells/data/"
io$hvfdir <- "~/datasets/scMET_ms/ecker2017/all_cells/hvf/hits/"
io$outdir <- "/Users/ckapoura/datasets/scMET_ms/ecker2017/all_cells/dimred_umap/data/"
if (!dir.exists(io$outdir)) { dir.create(io$outdir, recursive = TRUE) }
io$outprefix <- paste0(io$outdir, "/8k_k27ac_k4me1")


#############
## Options ##
#############
# Define genomic contexts
opts$anno <- list(
  "prom_2000_2000" = c("prom_2000_2000"),
  "distal_H3K27ac_cortex" = c("distal_H3K27ac_cortex"),
  "H3K4me1_cortex" = c("H3K4me1_cortex")
)
opts$min.counts <- 3  # Minimum # of CpGs per feature
opts$number.hvf <- 4000 # Number of highly variable features


##############################
## Load pre-computed parameter estimates
##############################
fit_dt <- opts$anno %>%
  map(~ fread(sprintf("%s/%s_vb.txt.gz", io$fitdir, .))) %>%
  rbindlist %>% .[,c("mu_median","gamma_median","epsilon_median",
                     "gauss_var", "binom_var", "Feature","anno")] %>%
  setnames(c("mu", "gamma", "epsilon", "gauss_var", "binom_var", "id", "anno"))

###########################
## Load methylation data ##
###########################
met_dt <- opts$anno %>%
  map(~ fread(sprintf("%s/%s.tsv.gz",io$data_parsed,.), showProgress = FALSE,
              colClasses = c("character","character","factor","numeric","integer")) %>%
        .[V1%in%sample_metadata$sample]
  ) %>% rbindlist %>% setnames(c("sample","id","anno","rate","Ntotal"))
# Filter features by number of CpGs
met_dt <- met_dt[Ntotal >= opts$min.counts]
# Filter features by overdispersion
hvfs <- lapply(X = opts$anno, function(an) {
  hvf_hits <- fread(sprintf("%s/hvf_%s_epsilon.txt.gz", io$hvfdir, an))
  hvfs <- hvf_hits %>% setorder(-tail_prob, -epsilon) %>% head(n = opts$number.hvf) %>% .$id
  return(hvfs)
})
met_dt <- met_dt[id %in% Reduce(c, hvfs)]

################
## Parse data ##
################
# Calculate M value from Beta value
met_dt[, m := log2(((rate/100) + 0.01) / (1 - (rate/100) + 0.01))]

# Prepare data for MOFA
met_dt <- met_dt %>% .[,c("sample","id","m")] %>% setnames(c("sample","feature","value"))

####################
## Fit MOFA model ##
####################
object <- create_mofa(met_dt)
data_opts <- get_default_data_options(object)
model_opts <- get_default_model_options(object)
model_opts$likelihoods[1] <- "gaussian"
model_opts$num_factors <- 15
model_opts$ard_weights <- FALSE
model_opts$spikeslab_weights <- FALSE
train_opts <- get_default_training_options(object)
train_opts$convergence_mode <- "fast"
object <- prepare_mofa(
  object = object,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
model <- run_mofa(object, sprintf("%s.hdf5", io$outprefix), save_data = FALSE)
#model <- load_model(sprintf("%s.hdf5", io$outprefix), load_data = FALSE)

##########################
## Extract MOFA factors ##
##########################
# Remove factors that explain very little variance
factors <- names(which(model@cache$variance_explained$r2_per_factor[[1]][,1] > 0.0001))
model <- subset_factors(model, factors)

# Fetch factors and merge with metadat
Z.mofa <- get_factors(model)[[1]]
Z.mofa <- Z.mofa %>% as.data.table %>% .[,sample := rownames(Z.mofa)]


################
## Clustering ##
################
set.seed(12345)
###############################
## Get clustering statistics ##
###############################
cluster_purity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

# Perform multiple clusterings and take average performance
ari = purity <- c()
for (i in 1:300) {
  # Merge
  dt <- sample_metadata %>% merge(Z.mofa, by = "sample")
  # Define true labels
  true_clusters <- as.factor(dt$Neuron_type3)
  ntrue_clusters <- length(unique(true_clusters))
  # Run k-means clustering on the latent space defined by MOFA
  clustering <- cluster_samples(model, k = ntrue_clusters,
                                iter.max = 50)$cluster %>% as.factor
  dt <- dt %>% merge(data.table(sample = names(clustering),
                                cluster = clustering), by = "sample")
  ari[i] <- adjustedRandIndex(x = as.numeric(true_clusters),
                              y = as.numeric(dt$cluster))
  purity[i] <- cluster_purity(true_clusters, dt$cluster) %>% round(3)
}

# Store results
cl_res <- data.table(
  anno = paste(opts$anno, collapse = " "),
  model = opts$model,
  hvf = opts$number.hvf,
  ari = mean(ari),
  purity = mean(purity)
)

#################
## Save output ##
#################
fwrite(cl_res, sprintf("%s_clustering.txt", io$outprefix),
       sep = "\t", quote = TRUE)

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

# Source settings
source("../../load_settings.R")

io$fitdir <- "~/datasets/scMET_ms/ecker2017/all_cells/data/"
io$hvfdir <- "~/datasets/scMET_ms/ecker2017/all_cells/hvf/hits/"
opts$min.counts <- 3  # Minimum # of CpGs per feature

######################
## Define arguments ##
######################
p <- ArgumentParser(description = '')
p$add_argument('--model',     type = "character", help = 'Model')
p$add_argument('--anno',      type = "character", help = 'Genomic context')
p$add_argument('--hvf',       type = "integer",   help = 'Number of HVFs')
p$add_argument('--outprefix', type = "character", help = 'Output directory')
args <- p$parse_args(commandArgs(TRUE))

##############################
## Load pre-computed parameter estimates
##############################
fit_dt <- args$anno %>%
  map(~ fread(sprintf("%s/%s_vb.txt.gz", io$fitdir, .))) %>%
  rbindlist %>% .[,c("mu_median","gamma_median","epsilon_median",
                     "gauss_var", "binom_var", "Feature","anno")] %>%
  setnames(c("mu", "gamma", "epsilon", "gauss_var", "binom_var", "id", "anno"))

###########################
## Load methylation data ##
###########################
met_dt <- args$anno %>%
  map(~ fread(sprintf("%s/%s.tsv.gz",io$data_parsed,.), showProgress = FALSE,
              colClasses = c("character","character","factor","numeric","integer")) %>%
  .[V1%in%sample_metadata$sample]
) %>% rbindlist %>% setnames(c("sample","id","anno","rate","Ntotal"))
# Filter features by number of CpGs
met_dt <- met_dt[Ntotal >= opts$min.counts]

# Filter features by variability
if (args$model == "binomial") {
  hvfs <- fit_dt %>% setorder(-binom_var) %>% head(n = args$hvf) %>% .$id
} else if (args$model == "gaussian") {
  hvfs <- fit_dt %>% setorder(-gauss_var) %>% head(n = args$hvf) %>% .$id
} else if (args$model == "scmet") {
  hvf_hits <- fread(sprintf("%s/hvf_%s_epsilon.txt.gz", io$hvfdir, args$anno))
  hvfs <- hvf_hits %>% setorder(-tail_prob, -epsilon) %>% head(n = args$hvf) %>% .$id
} else {
  stop("Wrong model specification")
}
met_dt <- met_dt[id %in% hvfs]

################
## Parse data ##
################
# Calculate M value from Beta value
met_dt[, m := log2(((rate/100) + 0.01) / (1 - (rate/100) + 0.01))]

# Prepare data for MOFA
met_dt <- met_dt %>% .[,c("sample","id","m")] %>%
  setnames(c("sample","feature","value"))

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
model <- run_mofa(object, sprintf("%s.hdf5", args$outprefix), save_data = FALSE)
#model <- load_model(sprintf("%s.hdf5", args$outprefix), load_data = FALSE)

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
  anno = paste(args$anno, collapse = " "),
  model = args$model,
  hvf = args$hvf,
  ari = mean(ari),
  purity = mean(purity)
)

#################
## Save output ##
#################
fwrite(cl_res, sprintf("%s_clustering.txt", args$outprefix),
       sep = "\t", quote = TRUE)

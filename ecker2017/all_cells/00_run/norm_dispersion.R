# Compute normalized dispersion values similar to scRNA-seq studies,
# e.g. ....
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
set.seed(123)

# --------------------------------------------------
# get variable genes from normalized UMI counts
# Adapted from
# https://github.com/10XGenomics/single-cell-3prime-paper/blob/master/pbmc68k_analysis/util.R
# --------------------------------------------------
# m: Cells x Genes matrix of methylation levels (either M or Beta values)
.get_norm_dispersion <-function(m) {
  df <- data.frame(mean = colMeans(m, na.rm = TRUE),
                   cv = apply(m, 2 ,sd, na.rm = TRUE) / colMeans(m, na.rm = TRUE),
                   var = apply(m, 2, var, na.rm = TRUE))
  df$dispersion <- with(df, var / mean)
  df$mean_bin <- with(df, cut(mean, breaks = c(-Inf, quantile(mean, seq(0.1, 1, 0.05)), Inf)))
  var_by_bin <- ddply(df, "mean_bin", function(x) {
    data.frame(bin_median = median(x$dispersion),
               bin_mad = mad(x$dispersion))
  })
  df$bin_disp_median <- var_by_bin$bin_median[match(df$mean_bin, var_by_bin$mean_bin)]
  df$bin_disp_mad <- var_by_bin$bin_mad[match(df$mean_bin, var_by_bin$mean_bin)]
  df$dispersion_norm <- with(df, abs(dispersion - bin_disp_median) / bin_disp_mad)
  df$Feature <- colnames(m)
  return(df)
}


################
## Load settings
source("../../load_settings.R")

# annos: "distal_H3K27ac_cortex", "H3K4me1_cortex", "prom_2000_2000"
anno <- "prom_2000_2000"
outdir <- "~/datasets/scMET_ms/ecker2017/all_cells/data/"
is_test <- FALSE

####################################
## Load / filter methylation data ##
cat("Loading data...\n")
Y <- read_filter_ecker_data(filename = sprintf("%s/%s.tsv.gz", io$data_parsed, anno),
                            opts = opts, is_differential = FALSE)
# Subset features in testing mode
if (is_test) {
  print("Test mode activated: subsetting features")
  Y <- Y[Feature %in% head(unique(Y$Feature), n = 200)]
}
# Calculate M value from Beta value (for now leave to Beta values)
Y <- Y %>%
  .[, rate := met_reads / total_reads] %>%
  #.[, m := log2(((rate) + 0.01) / (1 - (rate) + 0.01))] %>%
  .[, m := rate] %>%
  .[, c("Feature", "Cell", "m")]

Y_mat <- dcast(Y, Feature ~ Cell, value.var = "m")
# Convert to data.frame and make 1st column as rownames
Y_mat <- Y_mat %>% as.data.frame %>% remove_rownames %>%
  column_to_rownames(var = "Feature") %>% t

# Compute normalized dispersion estimates
df <- .get_norm_dispersion(Y_mat)
df$anno <- anno

##########
## Save ##
cat("Storing results ...\n")
if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
fwrite(df, sprintf("%s/%s_norm_dispersion_beta.txt.gz", outdir, anno), sep = "\t")
cat("Finished!!\n")

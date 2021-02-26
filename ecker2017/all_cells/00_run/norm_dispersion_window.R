# Compute normalized dispersion values similar to scRNA-seq studies,
# as in Zheng et al. (2017), Massively parallel digital transcriptional
#         profiling of single cells, Nature Communications.
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
  # Potentially a bug in the original implementation. Removing absolute value
  df$dispersion_norm_sign <- with(df, (dispersion - bin_disp_median) / bin_disp_mad)
  df$Feature <- colnames(m)
  return(df)
}


################
## Load settings
source("../../load_settings.R")

# annos: "distal_H3K27ac_cortex", "H3K4me1_cortex", "prom_2000_2000" "window20000_step20000"
anno <- "window20000_step20000"
outdir <- "~/datasets/scMET_ms/ecker2017/all_cells/data/"
is_test <- FALSE
io$basedir <- "~/datasets/ecker2017_local/mouse/"
io$data_parsed <- paste0(io$basedir,"/parsed_CG")

####################################
## Load / filter methylation data ##
cat("Loading data...\n")
Y <- readRDS(file = sprintf("%s/%s.rds", io$data_parsed, anno))
#Y <- read_filter_ecker_data(filename = sprintf("%s/%s.tsv.gz", io$data_parsed, anno),
#                            opts = opts, is_differential = FALSE)
met_groups <- c("Feature", "anno")

#############################
## Read methylation data   ##
#############################
Y <- Y %>% .[sample %in% opts$cells] %>%
  setnames(c("Cell", "Feature", "anno", "rate", "total_reads")) %>%
  .[, met_reads := round((rate/100) * total_reads)]
cat("Total # cells: ", length(unique(Y$Cell)), "\n")
cat("Total # features: ", length(unique(Y$Feature)), "\n")

#############################
## Filter methylation data ##
#############################
cat("Filtering data...\n")
cat("By number of CpGs and minimum number of cells...\n")
# Filter features by number of CpGs
Y <- Y[total_reads >= opts$min.cpgs]
# Filter features by minimum number of cells
Y <- Y[,Ncells := .N, by = met_groups] %>%
  .[Ncells >= opts$min.cells] %>% .[,Ncells := NULL]
cat("Total # cells: ", length(unique(Y$Cell)), "\n")
cat("Total # features: ", length(unique(Y$Feature)), "\n")

# Keep regions that are not in mean methylation extremes
cat("Removing highly/lowly mean methylated features...\n")
Y <- Y[, mean_rate := mean(rate/100), by = c("Feature", "anno")] %>%
  .[mean_rate > opts$met_rate_low & mean_rate < opts$met_rate_high] %>%
  .[, mean_rate := NULL]
cat("Total # cells: ", length(unique(Y$Cell)), "\n")
cat("Total # features: ", length(unique(Y$Feature)), "\n")

# Remove non-variable features
cat("Removing non-variable features below", opts$var.thresh,  "...\n")
Y <- Y[, var := var(rate/100), by = c("Feature", "anno")] %>%
  .[, .SD[var > opts$var.thresh], by = "anno"] %>% .[,var := NULL]
cat("Total # cells: ", length(unique(Y$Cell)), "\n")
cat("Total # features: ", length(unique(Y$Feature)), "\n")

## Remove again features by minimum number of cells
# Filter features by minimum number of cells
cat("By minimum number of cells (again) ...\n")
Y <- Y[,Ncells := .N, by = met_groups] %>%
  .[Ncells >= opts$min.cells] %>% .[,Ncells := NULL]
cat("Total # cells: ", length(unique(Y$Cell)), "\n")
cat("Total # features: ", length(unique(Y$Feature)), "\n")

############################
## Parse methylation data ##
############################
Y <- Y[, c("Feature", "Cell", "total_reads", "met_reads")]
# Reorder data by feature names.
Y <- Y[order(Feature), ]


# Subset features in testing mode
if (is_test) {
  print("Test mode activated: subsetting features")
  Y <- Y[Feature %in% head(unique(Y$Feature), n = 200)]
}
# Compute methylation rate
Y <- Y %>%
  .[, m := met_reads / total_reads] %>%
  .[, c("Feature", "Cell", "m")]

# Convert to data.frame and make 1st column as rownames
Y_mat <- dcast(Y, Feature ~ Cell, value.var = "m")
Y_mat <- Y_mat %>% as.data.frame %>% remove_rownames %>%
  column_to_rownames(var = "Feature") %>% t

# Compute normalized dispersion estimates as in Zheng et al. (2017)
df <- .get_norm_dispersion(Y_mat)
df$anno <- anno

##########
## Save ##
cat("Storing results ...\n")
if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
fwrite(df, sprintf("%s/%s_normdisp.txt.gz", outdir, anno), sep = "\t")
cat("Finished!!\n")


###########
## Tests ##

# tt <- df[rownames(df) %in% c("distal_H3K27ac_cortex_11725", "distal_H3K27ac_cortex_13166",
#                              "distal_H3K27ac_cortex_28573", "distal_H3K27ac_cortex_1480",
#                              "distal_H3K27ac_cortex_24013", "distal_H3K27ac_cortex_25912"), ]
# plot(df$mean, df$dispersion_norm)


# tmp_dt <- Y_mat[, colnames(Y_mat) %in% c("distal_H3K27ac_cortex_11725", "distal_H3K27ac_cortex_13166",
#                                          "distal_H3K27ac_cortex_28573", "distal_H3K27ac_cortex_1480",
#                                          "distal_H3K27ac_cortex_24013", "distal_H3K27ac_cortex_25912",
#                                          "distal_H3K27ac_cortex_10558")]
#
# # library(vioplot)
# # vioplot(tmp_dt, col="gold", srt = 45)
#
#
# Y_tmp <- copy(Y)
# Y_tmp <- Y_tmp %>% .[Feature %in% c("distal_H3K27ac_cortex_24869", "distal_H3K27ac_cortex_7927",
#                             "distal_H3K27ac_cortex_18656"), ]
# gg <- ggplot(Y_tmp, aes(x = Feature, y = m, fill = Feature)) +
#   geom_jitter(size = 0.8, alpha = 0.6, width = 0.25) +
#   geom_violin(alpha = 0.5) +
#   #scale_fill_manual(values = opts$colors3) +
#   labs(x = NULL, y = "Methylation", title = "Top three features called HVF by NormDisp") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(colour = "black", size = rel(1.3), hjust = 0.5),
#     axis.title.y = element_text(colour = "black", size = rel(1.1)),
#     axis.title.x = element_text(colour = "black", size = rel(1.1)),
#     axis.text.x = element_text(colour = "black", size = rel(0.7), angle=0, hjust = 0.5),
#     #axis.text.x = element_text(colour="black", size=rel(1.0), angle=30, hjust=1),
#     axis.text.y = element_text(colour = "black", size = rel(0.8)),
#     axis.text = element_text(colour = "black", size = rel(1.0)),
#     strip.text = element_text(colour = "black", size = rel(1.2)),
#     legend.position = "none"
#   )
#
# out_dir <- "~/datasets/scMET_ms/ecker2017/all_cells/hvf/comparison/hits/"
# if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
#
# pdf(paste0(out_dir, "top_hvf_normdisp.pdf"), width = 8, height = 3)
# print(gg)
# dev.off()
#
#
# Y_tmp <- copy(Y)
# Y_tmp <- Y_tmp %>% .[Feature %in% c("distal_H3K27ac_cortex_9891", "distal_H3K27ac_cortex_6245",
#                                     "distal_H3K27ac_cortex_24135", "distal_H3K27ac_cortex_6244"), ]
# gg <- ggplot(Y_tmp, aes(x = Feature, y = m, fill = Feature)) +
#   geom_jitter(size = 0.8, alpha = 0.6, width = 0.25) +
#   geom_violin(alpha = 0.5) +
#   #scale_fill_manual(values = opts$colors3) +
#   labs(x = NULL, y = "Methylation", title = "scMET HVF features not called by NormDisp") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(colour = "black", size = rel(1.3), hjust = 0.5),
#     axis.title.y = element_text(colour = "black", size = rel(1.1)),
#     axis.title.x = element_text(colour = "black", size = rel(1.1)),
#     axis.text.x = element_text(colour = "black", size = rel(0.7), angle=0, hjust = 0.5),
#     #axis.text.x = element_text(colour="black", size=rel(1.0), angle=30, hjust=1),
#     axis.text.y = element_text(colour = "black", size = rel(0.8)),
#     axis.text = element_text(colour = "black", size = rel(1.0)),
#     strip.text = element_text(colour = "black", size = rel(1.2)),
#     legend.position = "none"
#   )
#
# pdf(paste0(out_dir, "scmet_hvf_not_normdisp.pdf"), width = 8, height = 3)
# print(gg)
# dev.off()
#
#
# Y_tmp <- copy(Y)
# Y_tmp <- Y_tmp %>% .[Feature %in% c("distal_H3K27ac_cortex_15339", "distal_H3K27ac_cortex_9891",
#                                     "distal_H3K27ac_cortex_26439"), ]
# gg <- ggplot(Y_tmp, aes(x = Feature, y = m, fill = Feature)) +
#   geom_jitter(size = 0.8, alpha = 0.6, width = 0.25) +
#   geom_violin(alpha = 0.5) +
#   #scale_fill_manual(values = opts$colors3) +
#   labs(x = NULL, y = "Methylation", title = "Top three HVFs called by scMET") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(colour = "black", size = rel(1.3), hjust = 0.5),
#     axis.title.y = element_text(colour = "black", size = rel(1.1)),
#     axis.title.x = element_text(colour = "black", size = rel(1.1)),
#     axis.text.x = element_text(colour = "black", size = rel(0.7), angle=0, hjust = 0.5),
#     #axis.text.x = element_text(colour="black", size=rel(1.0), angle=30, hjust=1),
#     axis.text.y = element_text(colour = "black", size = rel(0.8)),
#     axis.text = element_text(colour = "black", size = rel(1.0)),
#     strip.text = element_text(colour = "black", size = rel(1.2)),
#     legend.position = "none"
#   )
#
# pdf(paste0(out_dir, "top_hvf_scmet.pdf"), width = 8, height = 3)
# print(gg)
# dev.off()

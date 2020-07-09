suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))

source("../../../load_settings.R")
######################
## Define arguments ##
######################
p <- ArgumentParser(description = '')
p$add_argument('--id',     type = "character", nargs = '+', help = 'feature id(s)')
p$add_argument('--anno',   type = "character",              help = 'genomic context')
p$add_argument('--outdir', type = "character",              help = 'Output directory')
args <- p$parse_args(commandArgs(TRUE))

if (is.null(args$outdir)) { stop("Define output directory.") }

# Define filtering criteria
opts$min.cpg <- 3

###############################
## Load DNA methylation data ##
###############################
data <- fread(sprintf("%s/%s.tsv.gz",io$data_parsed,args$anno), showProgress = FALSE) %>%
  setnames(c("sample","id","anno","rate","Ntotal")) %>% .[id %in% args$id]
# Filter by coverage
data <- data[Ntotal >= opts$min.cpg]

# Merge methylation data and sample metadata
data <- data %>% merge(sample_metadata, by = "sample")

data$Neuron_type3 <- factor(data$Neuron_type3,
                            levels = c("mDL-3", "mDL-2", "mDL-1", "mL6-2", "mL6-1", "mL5-2",
                                       "mL5-1", "mL4" , "mL2/3", "mNdnf-12", "mVip", "mSst-12", "mPv"))

# Read annotation data file
anno_dt <- fread(sprintf("%s/%s.bed.gz",io$features, args$anno)) %>%
  setnames(c("chr","start","end","strand","id","anno")) %>%
  .[,chr := paste0("chr",chr)] %>%
  .[id %in% args$id] %>%
  .[, locus := paste0(chr, ": ", start, " - ", end)]

###############
## Box plots ##
###############
for (i in args$id) {
  p1 <- ggplot(data[id == i], aes(x = Neuron_type1, y = rate/100, fill = Neuron_type1)) +
    geom_jitter(size = 0.8, alpha = 0.6, width = 0.25) +
    geom_violin(alpha = 0.5) +
    scale_fill_manual(values = opts$colors1) +
    labs(x = NULL, y = "Methylation", title = NULL) +
    # labs(x = NULL, y="Methylation", title =  anno_dt[id == i, locus]) +
    theme_classic() +
    theme(
      plot.title = element_text(colour = "black", size = rel(0.7), hjust = 0.5),
      axis.title.y = element_text(colour = "black", size = rel(1.1)),
      axis.title.x = element_text(colour = "black", size = rel(1.1)),
      axis.text.x = element_text(colour = "black", size = rel(1.0)),
      #axis.text.x = element_text(colour="black", size=rel(1.0), angle=30, hjust=1),
      axis.text.y = element_text(colour = "black", size = rel(0.8)),
      axis.text = element_text(colour = "black", size = rel(1.0)),
      strip.text = element_text(colour = "black", size = rel(1.2)),
      legend.position = "none"
    )
  pdf(sprintf("%s/boxplot_%s.pdf",args$outdir,i),
      width = 2.7, height = 2, useDingbats = FALSE)
  print(p1)
  dev.off()

  p2 <- ggplot(data[id == i], aes(x = Neuron_type3, y = rate/100, fill = Neuron_type3)) +
    geom_jitter(size = 1, alpha = 0.8, width = 0.25) +
    geom_violin(alpha = 0.5, width = 1.0) +
    facet_grid(~Neuron_type1, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = opts$colors3) +
    #labs(x = NULL, y = "Methylation", title = NULL) +
    labs(x = NULL, y = "Methylation", title =  anno_dt[id == i, locus]) +
    theme_classic() +
    theme(
      plot.title = element_text(colour = "black", size = rel(0.8), hjust = 0.5),
      axis.title.y = element_text(colour = "black", size = rel(1.1)),
      axis.text.x = element_text(colour = "black", size = rel(0.9)),
      #axis.text.x = element_text(colour="black", size=rel(1.0), angle=30, hjust=1),
      axis.text.y = element_text(colour = "black", size = rel(0.9)),
      strip.text = element_text(colour = "black", size = rel(1.1)),
      legend.position = "none"
    )
  pdf(sprintf("%s/boxplot_%s_subpop.pdf",args$outdir,i),
      width = 7.6, height = 2.4, useDingbats = FALSE)
  print(p2)
  dev.off()
}

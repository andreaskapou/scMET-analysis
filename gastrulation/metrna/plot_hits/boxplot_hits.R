suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(BASiCS))
suppressPackageStartupMessages(library(ggplot2))

source("../../load_settings.R")

######################
## Define arguments ##
######################
p <- ArgumentParser(description = '')
p$add_argument('--gene',          type = "character",              help = 'Feature name for RNA expression (in ENSEMBL ID)')
p$add_argument('--met.id',        type = "character",              help = 'Feature name for DNA methylation')
p$add_argument('--met.anno',      type = "character",              help = 'Genomic context for DNA methylation')
p$add_argument('--stage_lineage', type = "character", nargs = '+', help = 'stage and lineages to plot')
p$add_argument('--datadir',       type = "character",              help = 'Data directory for BASiCS output')
p$add_argument('--outdir',        type = "character",              help = 'Output directory')
args <- p$parse_args(commandArgs(TRUE))

if (is.null(args$outdir)) { stop("Define output directory.") }

####################
## Define options ##
####################

# Define stages to plot
# Define colors for the omics
opts$color <- c(
  "RNA expression" = "gray80",
  "DNA methylation" = "gray80"
)

# Define cells to use
opts$met_cells <- sample_metadata %>% .[stage_lineage %in% args$stage_lineage, id_met]
opts$rna_cells <- sample_metadata %>% .[stage_lineage %in% args$stage_lineage, id_rna]

###############
## Load data ##
###############

# Load DNA methylation data
met <- fread(sprintf("%s/%s.tsv.gz", io$met_data_parsed, args$met.anno)) %>%
  setnames(c("id_met","id","anno","Nmet","N","value")) %>%
  .[, value := value / 100] %>%
  .[id %in% args$met.id] %>% .[N >= opts$min.cpgs]

# Load RNA data
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
# Load BASiCS object
basics_obj <- readRDS(file = paste0(args$datadir, "/rna_basics.rds"))
denoised <- BASiCS_DenoisedCounts(gex, basics_obj)
denoised <- log(denoised + 1)

# Gene metadata
meta <- data.table(ens_id = rownames(sce), rna_id = rowData(sce)$symbol)
rna <- denoised[which(rownames(denoised) ==
                       meta[which(meta$rna_id == args$gene), ]$ens_id), , drop = FALSE] %>%
  t %>% as.data.table(keep.rownames = "id_rna") %>%
  melt(id.vars = "id_rna", value.name = "value", variable.name = "id")
rna$id <- args$gene


# sce <- readRDS(io$rna)
# sce <- sce[rowData(sce)$symbol == args$gene, ]
# rna <- denoised %>% t %>% as.data.table(keep.rownames = "id_rna") %>%
#   melt(id.vars = "id_rna", value.name = "value", variable.name = "id")
# rna$id <- args$gene

# Update sample metadata
sample_metadata <- sample_metadata %>%
  .[,c("sample", "id_met", "id_rna", "stage", "stage_lineage", "lineage10x_2")] %>%
  .[id_met %in% opts$met_cells | id_rna %in% opts$rna_cells]

# Merge data with sample metadata
met <- merge(met, sample_metadata, by = "id_met")
rna <- merge(rna, sample_metadata, by = "id_rna")

# bind in a single data table
to.plot <- do.call("rbind",list(
  rna[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>%
    .[,assay := "RNA expression"],
  met[,c("sample","id","stage","stage_lineage","lineage10x_2","value")] %>%
    .[,assay := "DNA methylation"]
)) %>% .[,stage_lineage := gsub("_"," ",stage_lineage)] %>%
  .[,assay := factor(assay, levels = c("RNA expression","DNA methylation"))]

##############
## Boxplots ##
##############
p <- ggplot(to.plot, aes(x = stage, y = value)) +
  facet_wrap(~assay, ncol = 1, scales = "free_y") +
  geom_jitter(size = 0.5, color = "black", width = 0.1, alpha = 0.8) +
  geom_violin(aes(fill = assay), alpha = 0.75, size = 0.4) +
  #geom_boxplot(aes(fill = assay), alpha = 0.6, outlier.shape = NA, width = 0.3, size = 0.25) +
  scale_fill_manual(values = opts$color) +
  scale_color_manual(values = "black") +
  scale_x_discrete(expand = c(0, -1)) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_classic() +
  theme(
    axis.title.y = element_text(colour = "black", size = rel(1.1), vjust = 1.5),
    axis.text.x = element_text(size = rel(1.4), color = "black"),
    axis.text.y = element_text(colour = "black",size = rel(1)),
    axis.line = element_line(colour = "black", size = rel(0.7)),
    axis.ticks.y = element_line(size = rel(0.7)),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

##########
## Save ##
##########
outfile <- sprintf("%s/boxplot_rna%s_met%s.pdf", args$outdir, args$gene, args$met.id)
pdf(outfile, useDingbats = FALSE, width = 3, height = 4.8)
print(p)
dev.off()

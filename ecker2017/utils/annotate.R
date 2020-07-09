###################################################################
##  Script to overlap bismark output files with genomic features ##
###################################################################

# This script overlaps the output bismark files (individual CpG sites)
# with genomic features such as promoters, gene bodies, etc.

# - Preprocessing of annotations: collect all CpG sites
# - Preprocessing of samples: collect all CpG sites from mm10 using the package "BSgenome.Mmusculus.UCSC.mm10"
# - Annotate samples with the preprocessed annotations

## Input ##
# output bismark file with the following columns:
#  "chr","pos","strand","context","met_counts","total_counts","rate"

## Output ##
# (1) a tmp folder with a bunch of tsv files with the preprocessed samples and annotations.
#   For example:
#     tmp/[SAMPLE]_[FEATURE].tsv
# (2) H3K27ac.tsv:
#   sample	id	anno	rate	total_number_cpgs
# nuclei-261_S1_L001      CGI_10  CGI     0       15
# nuclei-261_S1_L001      CGI_100 CGI     0       22
# nuclei-261_S1_L001      CGI_10003       CGI     0       7
# nuclei-261_S1_L001      CGI_10018       CGI     0       8

suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
# suppressMessages(library(doParallel))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################
p <- ArgumentParser(description = '')
p$add_argument('-n','--cores', type = "integer" ,help = 'Number of cores', default = 1)

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

################
## Define I/O ##
################
source("../load_settings.R")

####################
## Define Options ##
####################
# Genomic annotations
# opts$annos <- "all"
opts$annos <- c("prom_2000_2000")
if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$features, pattern = "\\.bed$"),"\\.bed"),"[[", 1)
cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for genomic annotations: %s\n",io$features))
cat(sprintf("- Input folder for methylation files: %s\n",io$data_raw))
cat(sprintf("- Output folder: %s\n", io$data_parsed))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse = " ")))
# cat(sprintf("- Number of cores: %d\n",args$cores))
cat("\n")

##############################
## Load genomic annotations ##
##############################

anno_list <- opts$annos %>% map(~
  fread(
    sprintf("%s/%s.bed.gz",io$features,.),
    sep = "\t", header = FALSE, select = c(1,2,3,4,5), verbose = FALSE) %>%
  setnames(c("chr","start","end","strand","id")) %>%
  setkey(chr,start,end)
)
names(anno_list) <- opts$anno


################################################################
## Calculate DNA methylation rate per cell and genomic region ##
################################################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$data_parsed), recursive = TRUEs)

# Run in parallel
# registerDoParallel(cores=args$cores)
# invisible(foreach(i=1:length(sample_metadata$sample)) %dopar% {
for (i in 1:length(sample_metadata$sample)) {
  sample <- sample_metadata$sample[i]
  samples_processed <- list.files(sprintf("%s/tmp", io$data_parsed))
  if (sum(str_detect(samples_processed,str_c(sample,"_",opts$anno))) == length(opts$anno)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample))
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))
    # Read DNA methylation data
    dat_sample <- fread(sprintf("%s/%s.tsv.gz",io$data_raw,sample),
                        sep = "\t", verbose = FALSE, showProgress = FALSE) %>%
      setnames(c("chr","pos","strand","context","met_counts","total_counts","rate"))

    # Add 'start' and 'end' columns to do the overlap
    dat_sample <- dat_sample[,c("start","end") := list(pos,pos)][,pos := NULL] %>%
      .[,chr := as.factor(chr)] %>% setkey(chr,start,end)

    # Overlap data with genomic annotations
    for (anno in opts$anno) {
      fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$data_parsed,sample,anno)
      if (file.exists(paste0(fname.out,".gz"))) {
        cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))

        # Calculate methylation rate for each genomic region by summarising over all overlapping CG sites
        ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch = 0) %>%
          .[,"i.end" := NULL] %>% setnames("i.start", "pos") %>%
          .[,c("sample","anno") := list(sample,anno)] %>%
          .[,.(rate = round(mean(rate) * 100), weight = .N), keyby = .(sample, id, anno)]

        # Store and save results
        fwrite(ov, fname.out, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
        system(sprintf("gzip -f %s",fname.out))
      }
    }
  }
}
rm(anno_list)

# Concatenate and save
for (i in opts$anno) {
  tmp <- opts$cells %>% map(~ fread(sprintf("%s/tmp/%s_%s.gz",io$met_data_parsed,.,i))) %>% rbindlist
  fwrite(tmp, sprintf("%s/%s.tsv.gz", io$met_data_parsed,i), quote = FALSE, sep = "\t", col.names = FALSE)
}

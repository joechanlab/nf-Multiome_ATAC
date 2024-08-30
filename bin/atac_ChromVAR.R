#!/usr/bin/env Rscript

suppressMessages({
    library(ArchR)
    library(stringr)
    library(dplyr)
    library(reticulate)
    library(parallel)
    library(optparse)
})

# Configure
addArchRThreads(threads = 80)
addArchRGenome('hg38')

# Set up command line options
option_list <- list(
    make_option(c("-a", "--atac_dir"), type="character", default=NULL,
                help="Input directory with atac-seq results", metavar="CHARACTER"),
    make_option(c("-i", "--insilicochip_pkl"), type="character", default=NULL,
                help="Pickle file with in silico chip results", metavar="CHARACTER"),
    make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="Output directory", metavar="CHARACTER"),
    make_option(c("-p", "--python_path"), type="character",default=NULL,
                help="Path to Python executable", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$insilicochip_pkl) || is.null(opt$atac_dir) || is.null(opt$out_dir) || is.null(opt$python_path)) {
    print_help(opt_parser)
    stop("All arguments are required", call.=FALSE)
}

# Set Python path
set.seed(1)

# Load ArchR project
proj <- loadArchRProject(path = opt$atac_dir)

# Load initial TF scores
Sys.setenv(RETICULATE_PYTHON = opt$python_path)
pd <- import("pandas")
tf_peaks_idx <- pd$read_pickle(opt$insilicochip_pkl)

# Calculate ChromVAR
tf_peaks_df <- read.table(file.path(opt$atac_dir, 'peak_counts/peaks.csv'), sep=',', header=T)

tf_peaks <- list()
for (tf in names(tf_peaks_idx)) {
    ind <- tf_peaks_idx[[tf]] + 1
    select_df <- tf_peaks_df[ind,]
    tf_peaks[[tf]] <- GRanges(select_df$seqnames, IRanges(select_df$start, select_df$end))
}

proj <- addPeakAnnotations(proj, tf_peaks, name="TFs", force=TRUE)

proj <- addDeviationsMatrix(
    ArchRProj = proj,
    peakAnnotation = "TFs",
    force = TRUE
)

TFscores_meta <- getMatrixFromProject(proj, 'TFsMatrix')
TFscores_df <- as(assay(TFscores_meta), 'dgTMatrix')

# Save results
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = opt$out_dir, load = FALSE)

write.table(as.matrix(TFscores_df), file.path(opt$out_dir, 'tf_score.ChromVAR.tsv'), quote = FALSE, sep ='\t')

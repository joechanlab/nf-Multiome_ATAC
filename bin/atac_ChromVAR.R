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
    make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="Output directory", metavar="CHARACTER"),
    make_option(c("-p", "--python_path"), type="character",
                default="/usersoftware/chanj3/postprocessing/bin/python",
                help="Path to Python executable", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$out_dir)) {
    print_help(opt_parser)
    stop("--out_dir argument must be provided", call.=FALSE)
}

# Set Python path
Sys.setenv(RETICULATE_PYTHON = opt$python_path)

set.seed(1)

# Load ArchR project
proj <- loadArchRProject(path = opt$out_dir)

# Load initial TF scores
pd <- import("pandas")
annot_fn <- file.path(opt$out_dir, 'Annotations/tf_peaks.ISChIP_idx.pkl')
tf_peaks_idx <- pd$read_pickle(annot_fn)

# Calculate ChromVAR
tf_peaks_df <- read.table(file.path(opt$out_dir, 'export/peak_counts/peaks.csv'), sep=',', row.names=1, header=T)

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
proj <- saveArchRProject(ArchRProj = proj)

write.table(as.matrix(TFscores_df), file.path(opt$out_dir, 'Annotations/tf_score.ChromVAR.tsv'), quote = FALSE, sep ='\t')

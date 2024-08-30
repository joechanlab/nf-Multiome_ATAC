#!/usr/bin/env Rscript

suppressMessages({
    library(ArchR)
    library(parallel)
    library(anndata)
    library(Rphenograph)
    library(optparse)
})

set.seed(1)

# Configure
addArchRThreads(threads = 12)
addArchRGenome('hg38')

get_overlaps <- function(peaks, motif_matches, tf) {
    # Input:
    #   peaks: GRanges object containing peak regions
    #   motif_matches: GRanges object containing motif match regions
    #   tf: character, name of the transcription factor

    # Assign input GRanges objects to variables
    gr1 <- motif_matches
    gr2 <- peaks

    # Find overlaps between peaks and motif matches
    hits <- findOverlaps(gr2, gr1, ignore.strand = TRUE)

    # Get indices of overlapping regions
    idx <- queryHits(hits)  # Indices of overlapping peaks
    idx2 <- subjectHits(hits)  # Indices of overlapping motif matches

    # Create a DataFrame with overlap information
    values <- unique(DataFrame(
        tf = tf,
        chr = seqnames(gr2)[idx],
        start = start(gr2)[idx],
        end = end(gr2)[idx],
        motifscore = gr1$score[idx2]
    ))

    # Convert DataFrame to data.frame
    values <- as.data.frame(values)

    return(values)
}

# Parse command line arguments
option_list <- list(
    make_option(c("-f", "--input_fragment"), type="character", default=NULL,
                help="Path to input fragment file", metavar="FILE"),
    make_option(c("-a", "--input_h5ad"), type="character", default=NULL,
                help="Path to input h5ad file", metavar="FILE"),
    make_option(c("-s", "--sample_name"), type="character", default=NULL,
                help="Sample name", metavar="STRING"),
    make_option(c("-o", "--output_dir"), type="character", default=".",
                help="Output directory", metavar="PATH")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input_fragment) || is.null(opt$input_h5ad) || is.null(opt$sample_name)) {
    print_help(opt_parser)
    stop("All arguments are required.", call.=FALSE)
}

# Set filtering parameters
filterTSS <- 9
filterFrags <- ifelse(grepl('RU581_LN', opt$sample_name), 5000, 2000)

# Create Arrow files
ArrowFiles <- createArrowFiles(
    inputFiles = opt$input_fragment,
    sampleNames = opt$sample_name,
    minTSS = filterTSS,
    minFrags = filterFrags,
    maxFrags = 1e+20,
    addTileMat = TRUE,
    addGeneScoreMat = FALSE,
    excludeChr = c('chrM')
)

# Create ArchR project
proj_name <- opt$output_dir
proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = paste0(proj_name, "_raw"),
    copyArrows = FALSE
)

# Load RNA data and subset cells
adata <- read_h5ad(opt$input_h5ad)
multiome_cells <- paste0(opt$sample_name, "#", rownames(adata))
multiome_cells <- intersect(multiome_cells, getCellNames(proj))
proj <- subsetArchRProject(proj, multiome_cells, proj_name, force=TRUE)

# Perform dimensionality reduction and clustering
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix",
    name = "IterativeLSI", scaleDims=FALSE, force=TRUE, varFeatures=25000)
var_features <- proj@reducedDims[['IterativeLSI']]$LSIFeatures

# Add gene score matrices
chrs <- getChromSizes(proj)
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrixFull', force=TRUE)

# Perform clustering and peak calling
lsi <- getReducedDims(proj, reducedDims = 'IterativeLSI')
res_phenograph <- Rphenograph(lsi, k=15)
clusters_phenograph <- membership(res_phenograph[[2]])
proj$Clusters <- paste0('C', as.character(clusters_phenograph))
proj <- addGroupCoverages(proj, groupBy='Clusters', maxFragmentLength=147)
proj <- addReproduciblePeakSet(proj, groupBy='Clusters')
proj <- addPeakMatrix(proj, maxFragmentLength=147, ceiling=10^9)

# Add motif annotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)
peaks <- getPeakSet(proj)
motif_matches <- readRDS(proj@peakAnnotation[['Motif']]$Positions)
df <- data.frame(matrix(ncol = 5, nrow = 0))
for (tf in names(motif_matches)) {
    if (length(motif_matches[[tf]]) > 0) {
        df <- rbind(df, get_overlaps(peaks, motif_matches[[tf]], tf))
    }
}
df$tf_name <- sapply(df$tf, function(x) strsplit(x, '_')[[1]][1])
df$peak_name <- paste0(df$chr, ':', df$start, '-', df$end)
df <- df %>%
    dplyr::arrange(desc(motifscore)) %>%
    dplyr::distinct(tf, peak_name, .keep_all = TRUE)
motif_dir <- dirname(proj@peakAnnotation[['Motif']]$Positions)
ofile <- file.path(proj_name, 'PeaksOverlapMotifs.csv')
write.table(df, ofile, row.names=FALSE, quote=FALSE, sep='\t')

# Extract peak counts
peak_counts_dir <- file.path(proj_name, "peak_counts")
dir.create(peak_counts_dir)
peak.counts <- getMatrixFromProject(proj, 'PeakMatrix')
chr_order <- sort(seqlevels(peaks))
reordered_features <- lapply(chr_order, function(chr) peaks[seqnames(peaks) == chr])
reordered_features <- Reduce("c", reordered_features)
names(reordered_features) <- sprintf("Peak%d", seq_along(reordered_features))
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, file.path(peak_counts_dir, "counts.mtx"))
write.csv(colnames(peak.counts), file.path(peak_counts_dir, "cells.csv"), quote=FALSE, row.names=FALSE)
write.csv(as.data.frame(reordered_features), file.path(peak_counts_dir, "peaks.csv"), quote=FALSE, row.names=FALSE)

# Save ArchR project
proj <- saveArchRProject(ArchRProj = proj)

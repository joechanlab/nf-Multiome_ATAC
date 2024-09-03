#!/usr/bin/env Rscript

suppressMessages({
    library(ArchR)
    library(parallel)
    library(stringr)
    library(Matrix)
    library(anndata)
    library(optparse)
})

set.seed(1)

# Configure
addArchRThreads(threads = 80)
addArchRGenome('hg38')

# Parse command line arguments
option_list <- list(
    make_option(c("-s", "--sample"), type="character", default=NULL,
                help="Sample name", metavar="CHARACTER"),
    make_option(c("-a", "--atac_dir"), type="character", default=NULL,
                help="Directory containing Arrow files", metavar="PATH"),
    make_option(c("-r", "--rna_h5"), type="character", default=NULL,
                help="Path to 10x feature matrix h5 file", metavar="FILE"),
    make_option(c("-d", "--rna_h5ad"), type="character", default=NULL,
                help="Path to 10x feature matrix h5 file", metavar="FILE"),
    make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="Output directory", metavar="PATH")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$atac_dir) || is.null(opt$rna_h5) || is.null(opt$rna_h5ad)) {
    print_help(opt_parser)
    stop("All arguments are required.", call.=FALSE)
}

proj <- loadArchRProject(path = opt$atac_dir)

# Load processed scRNAs-seq h5ad
rna_h5ad <- read_h5ad(opt$rna_h5ad)

# Check if cell names start with opt$sample followed by #
if (!all(startsWith(rna_h5ad$obs_names, paste0(opt$sample, "#")))) {
    # Prepend the prefix only to cell names that don't have it
    rna_h5ad$obs_names <- ifelse(
        startsWith(rna_h5ad$obs_names, paste0(opt$sample, "#")),
        rna_h5ad$obs_names,
        paste0(opt$sample, "#", rna_h5ad$obs_names)
    )
}

# Gene Activity
tmp <- import10xFeatureMatrix(
    input = c(opt$rna_h5),
    names = c(opt$sample)
)
feature_df <- rowData(tmp)
g2 <- intersect(rownames(feature_df), colnames(rna_h5ad))
feature_df <- feature_df[g2,]

tmp <- str_split(feature_df$interval,'[:-]')
tmp2 <- data.frame(do.call(rbind, tmp))
colnames(tmp2) <- c('chr','start','stop')
rowRanges <- makeGRangesFromDataFrame(tmp2)

# Create single cell experiment
seRNA <- SummarizedExperiment(list(counts=t(rna_h5ad[,g2]$X)), rowRanges=rowRanges)

inputFiles <- list.files(file.path(opt$atac_dir, 'ArrowFiles'), full.names=T)
names(inputFiles) <- gsub('.arrow','',basename(inputFiles))

# Add scRNA
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    useMatrix = 'GeneExpressionMatrix',
    reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.25,
    resolution = 1,
    returnLoops = FALSE
)

p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
p2g <- as.data.frame(p2g)

# Save
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
proj <- saveArchRProject(ArchRProj = proj, outputDirectory = opt$out_dir)

# Gene scores
write.csv(p2g, file.path(opt$out_dir, 'peak2genes.csv'), quote=FALSE)

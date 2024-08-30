#!/usr/bin/env Rscript

suppressMessages({
    library(ArchR)
    library(parallel)
    library(stringr)
    library(Matrix)
    library(optparse)
})

set.seed(1)

# Configure
addArchRThreads(threads = 80)
addArchRGenome('hg38')

# Parse command line arguments
option_list <- list(
    make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="Output directory", metavar="PATH"),
    make_option(c("-m", "--mtx_file"), type="character", default=NULL,
                help="Path to counts matrix file", metavar="FILE"),
    make_option(c("-f", "--feature_matrix"), type="character", default=NULL,
                help="Path to 10x feature matrix h5 file", metavar="FILE"),
    make_option(c("-a", "--arrow_dir"), type="character", default=NULL,
                help="Directory containing Arrow files", metavar="PATH")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$out_dir) || is.null(opt$mtx_file) || is.null(opt$feature_matrix) || is.null(opt$arrow_dir)) {
    print_help(opt_parser)
    stop("All arguments are required.", call.=FALSE)
}

setwd(opt$out_dir)

proj <- loadArchRProject(path = opt$out_dir)

# Gene Activity
counts <- t(as.matrix(read.table(opt$mtx_file, header=T, row.names=1, sep=',', comment.char = "!")))

tmp <- import10xFeatureMatrix(
    input = c(opt$feature_matrix),
    names = c("RU1042_PDX")
)

feature_df <- rowData(tmp)
g2 <- intersect(rownames(feature_df), rownames(counts))
feature_df <- feature_df[g2,]

tmp <- str_split(feature_df$interval,'[:-]')
tmp2 <- data.frame(do.call(rbind, tmp))
colnames(tmp2) <- c('chr','start','stop')
rowRanges <- makeGRangesFromDataFrame(tmp2)

# Create single cell experiment
seRNA <- SummarizedExperiment(list(counts=counts[g2,]), rowRanges=rowRanges)

inputFiles <- list.files(opt$arrow_dir, full.names=T)
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
proj <- saveArchRProject(ArchRProj = proj)

# Gene scores
write.csv(p2g, file.path(opt$out_dir, 'export/peak2genes.csv'), quote=FALSE)

#!/usr/bin/env Rscript

suppressMessages({
    library(Signac)
    library(Seurat)
    library(GenomeInfoDb)
    library(EnsDb.Hsapiens.v86)
    library(ggplot2)
    library(patchwork)
})

set.seed(1)

# Configure
addArchRThreads(threads = 12)
addArchRGenome('hg38')

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

# Input files
filterTSS = 9
filterFrags = ifelse(grepl('RU581_LN', opt$sample_name), 5000, 2000)

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

# Create project
proj_name <- opt$sample_name
proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = opt$output_dir,
    copyArrows = FALSE
)

# Load RNA files
# Subset of cells determined in RNA
adata = read_h5ad(opt$input_h5ad)
multiome_cells <- paste0(opt$sample_name, "#", rownames(adata))
multiome_cells <- intersect(multiome_cells, getCellNames(proj))

# Subset
proj <- subsetArchRProject(proj, multiome_cells, proj_name, force=T)

# SVD, Clustering, UMAP
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix",
    name = "IterativeLSI", scaleDims=FALSE, force=TRUE, varFeatures=25000)
var_features = proj@reducedDims[['IterativeLSI']]$LSIFeatures

# Gene scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)
# Full gene activity
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrixFull', force = T)

# Peaks
lsi = getReducedDims(proj, reducedDims = 'IterativeLSI')
res_phenograph = Rphenograph(lsi, k=15)
clusters_phenograph = membership(res_phenograph[[2]])
proj$Clusters = paste0('C', as.character(clusters_phenograph))
proj <- addGroupCoverages(proj, groupBy='Clusters', maxFragmentLength=147)
proj <- addReproduciblePeakSet(proj, groupBy='Clusters')
proj <- addPeakMatrix(proj, maxFragmentLength=147, ceiling=10^9)

# Add motif matches
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif_cisbp", force = T)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "vierstra", name = "Motif_vierstra", collection='archetype', force = T)

# Save
proj <- saveArchRProject(ArchRProj = proj)

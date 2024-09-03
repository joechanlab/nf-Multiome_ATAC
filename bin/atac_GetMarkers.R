#!/usr/bin/env Rscript

suppressMessages({
    library(gtools)
    library(HDF5Array)
    library(Seurat)
    library(Matrix)
    library(tidyr)
    library(dplyr)
    library(parallel)
    library(ArchR)
    library(anndata)
    library(optparse)
})

# Configure
addArchRThreads(threads = 20)
addArchRGenome('hg38')

# Parse command line arguments
option_list <- list(
    make_option(c("-s", "--sample"), type="character", default=NULL,
                help="Sample name", metavar="CHARACTER"),
    make_option(c("-a", "--atac_dir"), type="character", default=NULL,
                help="Directory containing Arrow files", metavar="PATH"),
    make_option(c("-d", "--rna_h5ad"), type="character", default=NULL,
                help="RNA h5ad file", metavar="PATH"),
    make_option(c("-g", "--group"), type="character", default="leiden",
                help="Variable to group by", metavar="CHARACTER"),
    make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="Output directory", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$sample) || is.null(opt$atac_dir) || is.null(opt$rna_h5ad) || is.null(opt$group) || is.null(opt$out_dir)) {
    print_help(opt_parser)
    stop("Both --sample and --atac_dir and --rna_h5ad and --group and --out_dir must be provided.", call.=FALSE)
}

# Load ArchR project
proj <- loadArchRProject(path = opt$atac_dir)
peaks <- getPeakSet(proj)

# Add the group to the sampleColData
rna_h5ad <- read_h5ad(opt$rna_h5ad)
group_labels <- rna_h5ad$obs[[opt$group]]
names(group_labels) <- rownames(rna_h5ad$obs)

# Get the cell names from the ArchR project
atac_cells <- getCellNames(proj)

# Match the cell names between RNA and ATAC data
# Assuming the cell names in ATAC data have the sample name as prefix
atac_cells_stripped <- sub("^[^#]+#", "", atac_cells)
matching_indices <- match(atac_cells_stripped, names(group_labels))

# Create a new vector with matched group labels for ATAC cells
atac_group_labels <- as.character(group_labels[matching_indices])
proj <- addCellColData(proj, data = atac_group_labels, name = 'group', cells = atac_cells, force = TRUE)

# Validate that the labels were added correctly
if (length(proj$group) != length(getCellNames(proj))) {
    stop("Error: The number of group labels does not match the number of cells in the ArchR project.")
}
ct = unique(proj$group)

# Function to process markers
process_markers <- function(matrix_name, output_dir, file_prefix, sample) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (length(ct) >= 2) {
        comparisons <- permutations(length(ct), 2, ct)
        # Filter comparisons to ensure each group has at least 20 cells
        valid_comparisons <- comparisons[apply(comparisons, 1, function(comp) {
            all(table(proj$group[proj$group %in% comp]) >= 20)
        }), , drop = FALSE]
        if (nrow(valid_comparisons) == 0) {
            warning("No valid comparisons found with at least 20 cells in each group.")
            return()
        }
        comparisons <- valid_comparisons
        for (i in 1:nrow(comparisons)) {
            # Print the number of cells in each group for this comparison
            group1_cells <- sum(proj$group == comparisons[i,1])
            group2_cells <- sum(proj$group == comparisons[i,2])
            markers <- getMarkerFeatures(
                ArchRProj = proj,
                useMatrix = matrix_name,
                groupBy = "group",
                bias = c("TSSEnrichment", "log10(nFrags)"),
                testMethod = "wilcoxon",
                useGroups = c(comparisons[i,2]),
                bgdGroups = c(comparisons[i,1]),
                threads = 20
            )
            comp2 <- gsub(paste0(sample, ': '), '', comparisons[i,2])
            comp1 <- gsub(paste0(sample, ': '), '', comparisons[i,1])
            out_dir2 <- file.path(output_dir, sprintf('%s.%s_%s', sample, comp2, comp1))
            dir.create(out_dir2, showWarnings = FALSE, recursive = TRUE)
            saveHDF5SummarizedExperiment(markers,
                                        dir = out_dir2,
                                        prefix = sprintf('%s.%s.%s_%s.', file_prefix, sample, comp2, comp1),
                                        replace = TRUE)

            markerList <- getMarkers(markers, cutOff = "FDR <= 1 & Log2FC > 0")
            out_file <- file.path(out_dir2, sprintf('%s.%s.%s_%s.csv', file_prefix, sample, comp2, comp1))
            write.table(markerList, out_file, sep = '\t', quote = FALSE)
        }
    }
}

# Process markers
process_markers("PeakMatrix", sprintf('%s/marker_peaks/%s/', opt$out_dir, opt$group), "MarkerPeaks", opt$sample)
process_markers("TFsMatrix", sprintf('%s/marker_TFs/%s/', opt$out_dir, opt$group), "MarkerTFs", opt$sample)
process_markers("GeneScoreMatrixFull", sprintf('%s/marker_genes/%s/', opt$out_dir, opt$group), "MarkerGenes", opt$sample)

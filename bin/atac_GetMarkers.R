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
    library(optparse)
})

# Configure
addArchRThreads(threads = 20)
addArchRGenome('hg38')

# Parse command line arguments
option_list <- list(
    make_option(c("-s", "--sample_query"), type="character", default=NULL,
                help="Sample query", metavar="CHARACTER"),
    make_option(c("-c", "--cell_type_category"), type="character", default=NULL,
                help="Cell type category", metavar="CHARACTER"),
    make_option(c("-o", "--out_dir"), type="character", default="/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.ATAC.metacells.combined.042023/ArchR/",
                help="Output directory", metavar="CHARACTER"),
    make_option(c("-r", "--rna_dir"), type="character", default="/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.individual.042023/",
                help="RNA directory", metavar="CHARACTER"),
    make_option(c("-t", "--target_mc"), type="integer", default=25,
                help="Target metacell number", metavar="INTEGER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$sample_query) || is.null(opt$cell_type_category)) {
    print_help(opt_parser)
    stop("Both --sample_query and --cell_type_category must be provided.", call.=FALSE)
}

# Set working directory
setwd(opt$out_dir)

# Load ArchR project
proj <- loadArchRProject(path = opt$out_dir)

peaks <- getPeakSet(proj)

sample_labels <- rownames(proj@sampleColData)

# Process cell types
cell_type <- list()
for (sample in sample_labels) {
    print(sample)
    fn <- file.path(opt$rna_dir, sample, sprintf('obs.%s.042023.csv', sample))
    df <- read.table(fn, sep = '\t', header=T, row.names=1)
    num_bc <- nrow(df) - 1
    num_mc <- floor(num_bc/opt$target_mc)
    if (sample %in% c('RU599_PDX', 'RU298_PDX')) {
        fn_check <- sprintf('%s%s/SEACells/%d/SEACell_obs.%s.%d.042023.csv', opt$rna_dir, sample, num_mc, sample, num_mc)
        df2 <- read.table(fn_check, sep='\t', header=T, row.names=1)
        rownames(df2) <- paste0(sample, '#', rownames(df2))
        cell_type <- rbind(cell_type, df2 %>% dplyr::select(`cell_type_category`) %>% dplyr::rename(cell_type=`cell_type_category`))
    } else {
        fn_check <- Sys.glob(sprintf('%s%s/SEACells/%d/DomainAdapt_obs.%s.*.csv', opt$rna_dir, sample, num_mc, sample))
        if (length(fn_check) != 0) {
            df2 <- read.table(fn_check, sep=',', header=T, row.names=1)
            df2 <- df2 %>% dplyr::filter(Modality=='Paired')
            rownames(df2) <- gsub('Paired_', paste0(sample, '#'), rownames(df2))
            cell_type <- rbind(cell_type, df2 %>% dplyr::select(`cell_type_category`) %>% dplyr::rename(cell_type=`cell_type_category`))
        }
    }
}

cell_type <- cell_type[proj$cellNames,]

proj$cell_type <- cell_type

proj$sample_celltype <- paste0(proj$Sample, ': ', proj$cell_type)

s_ct <- levels(factor(proj$sample_celltype))

ct <- s_ct[grepl(opt$sample_query, s_ct)]

# Function to process markers
process_markers <- function(matrix_name, output_dir, file_prefix) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    if (length(ct) >= 2) {
        comparisons <- permutations(length(ct), 2, ct)
        for (i in 1:nrow(comparisons)) {
            markers <- getMarkerFeatures(
                ArchRProj = proj,
                useMatrix = matrix_name,
                groupBy = "sample_celltype",
                bias = c("TSSEnrichment", "log10(nFrags)"),
                testMethod = "wilcoxon",
                useGroups = c(comparisons[i,2]),
                bgdGroups = c(comparisons[i,1]),
                threads = 20
            )

            comp2 <- gsub(paste0(opt$sample_query, ': '), '', comparisons[i,2])
            comp1 <- gsub(paste0(opt$sample_query, ': '), '', comparisons[i,1])

            out_dir2 <- file.path(output_dir, sprintf('%s.%s_%s', opt$sample_query, comp2, comp1))
            dir.create(out_dir2, showWarnings = FALSE, recursive = TRUE)

            saveHDF5SummarizedExperiment(markers,
                                        dir = out_dir2,
                                        prefix = sprintf('%s.%s.%s_%s.', file_prefix, opt$sample_query, comp2, comp1),
                                        replace = TRUE)

            markerList <- getMarkers(markers, cutOff = "FDR <= 1 & Log2FC > 0")

            out_file <- file.path(out_dir2, sprintf('%s.%s.%s_%s.csv', file_prefix, opt$sample_query, comp2, comp1))
            write.table(markerList, out_file, sep = '\t', quote = FALSE)
        }
    }
}

# Process different types of markers
process_markers("PeakMatrix", sprintf('%s/export/marker_peaks/%s/', opt$out_dir, opt$cell_type_category), "MarkerPeaks")
process_markers("TFsMatrix", sprintf('%s/export/marker_TFs/%s/', opt$out_dir, opt$cell_type_category), "MarkerTFs")
process_markers("GeneScoreMatrixFull", sprintf('%s/export/marker_genes/%s/', opt$out_dir, opt$cell_type_category), "MarkerGenes")
process_markers("TFsSpearmanMatrix", sprintf('%s/export/marker_TFs_spearman/%s/', opt$out_dir, opt$cell_type_category), "MarkerTFs_spearman")
process_markers("TFsIncludeRepressMatrix", sprintf('%s/export/marker_TFs_include_repress/%s/', opt$out_dir, opt$cell_type_category), "MarkerTFs_include_repress")
process_markers("TFsSpearmanIncludeRepressMatrix", sprintf('%s/export/marker_TFs_spearman_include_repress/%s/', opt$out_dir, opt$cell_type_category), "MarkerTFs_spearman_include_repress")

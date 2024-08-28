#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATAC_QC } from './modules/atac_qc'
include { ATAC_PREPROCESS } from './modules/atac_preprocess'

workflow {
    // Access the samplesheet
    sample_sheet = file(params.samplesheet)

    // Read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // Create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def (name, rna_h5ad_path, atac_h5ad_path, fragment_path, fragment_index_path) = row
        return tuple(name, file(rna_h5ad_path), file(atac_h5ad_path), file(fragment_path), file(fragment_index_path))
    }

    // Run ATAC QC
    ATAC_QC(ch_input)

    // Run ATAC preprocessing
    ATAC_PREPROCESS(ch_input)
}

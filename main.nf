#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ATAC_QC } from './modules/atac_qc'
include { ATAC_PREPROCESS } from './modules/atac_preprocess'
include { ATAC_INSILICOCHIP } from './modules/atac_insilicochip'
include { ATAC_CHROMVAR } from './modules/atac_chromvar'
include { ATAC_PEAK2GENE } from './modules/atac_peak2gene'
include { ATAC_GETMARKERS } from './modules/atac_getmarkers'
include { ATAC_PREPROCESS_SNAPATAC2 } from './modules/atac_preprocess_SnapATAC2'
include { ATAC_JOINT_EMBEDDING } from './modules/atac_joint_embedding'
include { ATAC_NETWORK_ANALYSIS } from './modules/atac_network_analysis'
include { ATAC_REPORT } from './modules/atac_report'

workflow {
    // Access the samplesheet
    sample_sheet = file(params.samplesheet)

    // Read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // Create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def (name, rna_h5, rna_h5ad, rna_seacells_h5ad, rna_seacells_dir, atac_h5ad, fragment_path, fragment_index_path) = row
        return tuple(name, rna_h5, rna_h5ad, rna_seacells_h5ad, rna_seacells_dir, atac_h5ad, fragment_path, fragment_index_path)
    }

    // Run SnapATAC2 instead for preprocessing
    if (params.snapatac2) {
        // Run SnapATAC2 preprocessing
        ATAC_PREPROCESS_SNAPATAC2(ch_input)

        // Run joint embedding
        ATAC_JOINT_EMBEDDING(ATAC_PREPROCESS_SNAPATAC2.out.name, ATAC_PREPROCESS_SNAPATAC2.out.output_dir, ATAC_PREPROCESS_SNAPATAC2.out.rna_h5ad)

        // Run network analysis
        ATAC_NETWORK_ANALYSIS(ATAC_JOINT_EMBEDDING.out.name, ATAC_JOINT_EMBEDDING.out.output_dir)

        // Run report
        ATAC_REPORT(ATAC_NETWORK_ANALYSIS.out.name, ATAC_NETWORK_ANALYSIS.out.output_dir)

    } else {
        // Run ATAC QC
        ATAC_QC(ch_input)

        // Run ATAC preprocessing
        ATAC_PREPROCESS(ch_input)

        // Run In Silico CHIP
        ATAC_INSILICOCHIP(ATAC_PREPROCESS.out.name, ATAC_PREPROCESS.out.output_dir, ATAC_PREPROCESS.out.rna_seacells_h5ad, ATAC_PREPROCESS.out.rna_seacells_dir, ATAC_PREPROCESS.out.rna_h5, ATAC_PREPROCESS.out.rna_h5ad)

        // Run ChromVAR
        ATAC_CHROMVAR(ATAC_INSILICOCHIP.out.name, ATAC_INSILICOCHIP.out.arrow_dir, ATAC_INSILICOCHIP.out.output_pkl, ATAC_INSILICOCHIP.out.rna_h5, ATAC_INSILICOCHIP.out.rna_h5ad)

        // Run Peak2Gene
        ATAC_PEAK2GENE(ATAC_CHROMVAR.out.name, ATAC_CHROMVAR.out.out_dir, ATAC_CHROMVAR.out.rna_h5, ATAC_CHROMVAR.out.rna_h5ad)

        // ATAC_GETMARKERS
        ATAC_GETMARKERS(ATAC_PEAK2GENE.out.name, ATAC_PEAK2GENE.out.out_dir, ATAC_PEAK2GENE.out.rna_h5ad)
    }
}

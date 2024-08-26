#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {ATAC_QC} from './modules/atac_qc'

workflow {
    // access the samplesheet
    sample_sheet = file(params.samplesheet)

    // read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def name = row[0]
        def h5_path = file(row[1])
        def fragment_path = file(row[2])
        def fragment_index_path = file(row[3])
        return tuple(name, h5_path, fragment_path, fragment_index_path)
    }

    // run Cellbender
    ATAC_QC(ch_input)

}

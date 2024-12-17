#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_edta_pipeline'
include { EDTA                      } from './workflows/edta.nf'

workflow OUSHUJUN_EDTA {

    take:
    ch_genome

    main:
    EDTA ( ch_genome )
}

workflow {

    main:
    PIPELINE_INITIALISATION (
        params.version,
        args,
        params.outdir,
        params.genome
    )

    OUSHUJUN_EDTA ( PIPELINE_INITIALISATION.out.genome )
}
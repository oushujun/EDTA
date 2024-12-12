#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.genome           = 'genomes/*'
params.species          = 'others'
params.outdir           = 'results'

include { EDTA          } from './workflows/edta.nf'

workflow {
    EDTA()
}

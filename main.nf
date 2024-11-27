#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.genomes          = 'genomes/*'
params.species          = 'others'
params.cds              = ''
params.curatedlib       = ''
params.rmlib            = ''
params.sensitive        = false
params.anno             = false
params.rmout            = ''
params.maxdiv           = 40
params.evaluate         = true
params.exclude          = ''
params.maxint           = 5000
params.outdir           = 'results'

include { EDTA          } from './workflows/edta.nf'

// Test run: 
// ./main.nf -profile docker,test
// ./main.nf -profile conda,test

workflow {
    EDTA()
}

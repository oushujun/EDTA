#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { EDTA          } from './workflows/edta.nf'

workflow {
    EDTA()
}

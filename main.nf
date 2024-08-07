#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.genomes          = 'genomes/*' // To allow for more flexibility when specifying params.
// I'll be testing the whole pipeline with a single public chromosome from nf-core test datasets
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

// TODO: Check inputed repeat libraries, CDS, etc...
// TODO: Check exclude file

include { SANITIZE_HEADERS } from './modules/local/sanitize/main.nf'

// Test run: 
// ./main.nf --genomes https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta -profile docker
workflow {
    // - All nf-core pipelines/modules the [ [`meta`], data] pattern. I think we should follow that
    // so that our local and nf-core modules are interoperable
    //
    // - I am also adding a bit of personal code style, please rever it if you don't like it

    ch_genome                           = Channel.fromPath(params.genomes) // The channel emits a single genome at a time. That's why ch_genome

    // Create a meta object for each genome
    ch_meta_genome                         = ch_genome.map { genome -> 
        meta                            = [:]
        meta.id                         = genome.baseName
        return tuple(meta, genome)
    }

    // MODULE: sanitize - All modules and workflow names should follow CAPITAL_SNAKE

    SANITIZE_HEADERS ( ch_meta_genome )
}
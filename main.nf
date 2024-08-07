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

// Max resource options
params.max_cpus         = 12
params.max_memory       = '16.GB'
params.max_time         = '1.hour'

// TODO: Check inputed repeat libraries, CDS, etc...
// TODO: Check exclude file

include { SANITIZE_HEADERS } from './modules/local/sanitize/main.nf'

// Test run: 
// ./main.nf --genomes https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta -profile docker
// ./main.nf --genomes https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta -profile conda
workflow {

    ch_genome                           = Channel.fromPath(params.genomes)

    // Create a meta object for each genome
    ch_meta_genome                      = ch_genome.map { genome -> 
                                            meta        = [:]
                                            meta.id     = genome.baseName
                                            
                                            [ meta, genome ]
                                        }

    // MODULE: SANITIZE_HEADERS
    SANITIZE_HEADERS ( ch_meta_genome )
}
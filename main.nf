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

include { SANITIZE_HEADERS  } from './modules/local/sanitize/main.nf'
include { LTRHARVEST        } from './modules/nf-core/ltrharvest/main.nf'

// Test run: 
// ./main.nf -profile docker,test
// ./main.nf -profile conda,test
workflow {

    // Versions channel
    ch_versions                         = Channel.empty()

    
    ch_genome                           = Channel.fromPath(params.genomes)

    // Create a meta object for each genome
    ch_meta_genome                      = ch_genome.map { genome -> 
                                            meta        = [:]
                                            meta.id     = genome.baseName
                                            
                                            [ meta, genome ]
                                        }

    // MODULE: SANITIZE_HEADERS
    SANITIZE_HEADERS ( ch_meta_genome )

    ch_sanitized_fasta                  = SANITIZE_HEADERS.out.fasta

    // MODULE: LTRHARVEST
    LTRHARVEST ( ch_sanitized_fasta )

    ch_ltrharvest_gff3                  = LTRHARVEST.out.gff3
    ch_ltrharvest_scn                   = LTRHARVEST.out.scn

    ch_versions                         = ch_versions.mix(LTRHARVEST.out.versions)
}
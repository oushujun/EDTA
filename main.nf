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

include { SANITIZE_HEADERS          } from './modules/local/sanitize/main.nf'
include { LTRHARVEST                } from './modules/nf-core/ltrharvest/main.nf'
include { LTRFINDER                 } from './modules/nf-core/ltrfinder/main.nf'
include { CAT_CAT                   } from './modules/nf-core/cat/cat/main.nf'
include { LTRRETRIEVER_LTRRETRIEVER } from './modules/nf-core/ltrretriever/ltrretriever/main.nf'
include { TIRLEARNER                } from './modules/gallvp/tirlearner/main.nf'
// nf-core -v modules -g https://github.com/GallVp/nxf-components.git install
include { ANNOSINE                  } from './modules/gallvp/annosine/main.nf'

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

    // MODULE: LTRFINDER
    LTRFINDER  { ch_sanitized_fasta }

    ch_ltrfinder_gff3                   = LTRFINDER.out.gff
    ch_ltrfinder_scn                    = LTRFINDER.out.scn

    ch_versions                         = ch_versions.mix(LTRFINDER.out.versions)

    // MODULE: CAT_CAT
    ch_cat_cat_inputs                   = ch_ltrharvest_scn
                                        | join(ch_ltrfinder_scn)
                                        | map { meta, harvested, found -> [ meta, [ harvested, found ] ] }
    CAT_CAT ( ch_cat_cat_inputs )

    ch_ltr_candidates                   = CAT_CAT.out.file_out
    ch_versions                         = ch_versions.mix(CAT_CAT.out.versions.first())

    // MODULE: LTRRETRIEVER_LTRRETRIEVER
    ch_ltrretriever_inputs              = ch_sanitized_fasta.join(ch_ltr_candidates)

    LTRRETRIEVER_LTRRETRIEVER (
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> [ meta, fasta ] },
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> ltr },
        [],
        [],
        []
    )

    ch_ltrretriever_log             = LTRRETRIEVER_LTRRETRIEVER.out.log
    ch_pass_list                    = LTRRETRIEVER_LTRRETRIEVER.out.pass_list
    ch_annotation_out               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_out
    ch_annotation_gff               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_gff
    ch_ltrlib                       = LTRRETRIEVER_LTRRETRIEVER.out.ltrlib
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LTRRETRIEVER.out.versions.first())

    // MODULE: TIRLEARNER
    TIRLEARNER (
        ch_sanitized_fasta,
        params.species
    )

    ch_tirlearner_filtered_gff          = TIRLEARNER.out.filtered_gff
    ch_versions                         = ch_versions.mix(TIRLEARNER.out.versions)

    // These can also run in parallel
    // MODULE: ANNOSINE
    ANNOSINE (
        ch_sanitized_fasta,
        3 // mode
    )

    // Currently it's a topic, so need to fix that
    ch_versions                         = ch_versions.mix(ANNOSINE.out.versions)
    cb_annosine_seed_sine               = ANNOSINE.out.fa

}

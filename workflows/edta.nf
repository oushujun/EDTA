include { SANITIZE_HEADERS                  } from '../modules/local/sanitize/main.nf'
include { LTRHARVEST                        } from '../modules/nf-core/ltrharvest/main.nf'
include { LTRFINDER                         } from '../modules/nf-core/ltrfinder/main.nf'
include { CAT_CAT                           } from '../modules/nf-core/cat/cat/main.nf'
include { LTRRETRIEVER_LTRRETRIEVER         } from '../modules/nf-core/ltrretriever/ltrretriever/main.nf'
include { TIRLEARNER                        } from '../modules/gallvp/tirlearner/main.nf'
include { ANNOSINE                          } from '../modules/gallvp/annosine/main.nf'
include { REPEATMODELER_BUILDDATABASE       } from '../modules/nf-core/repeatmodeler/builddatabase/main.nf'
include { REPEATMODELER_REPEATMODELER       } from '../modules/nf-core/repeatmodeler/repeatmodeler/main.nf'
include { FASTA_HELITRONSCANNER_SCAN_DRAW   } from '../subworkflows/gallvp/fasta_helitronscanner_scan_draw/main.nf'
include { FORMAT_HELITRONSCANNER_OUT        } from '../modules/local/format_helitronscanner_out/main.nf'
include { LTR_RETRIEVER_POSTPROCESS         } from '../modules/local/ltr_retriever_postprocess/main.nf'
include { TIR_LEARNER_POSTPROCESS           } from '../modules/local/tir_learner_postprocess/main.nf'

include { softwareVersionsToYAML            } from '../modules/local/utils/main.nf'

workflow EDTA {

    main:

    // Versions channel
    ch_versions                         = Channel.empty()

    
    ch_genome                           = Channel.fromPath(params.genomes)

    // Create a meta object for each genome
    ch_meta_genome                      = ch_genome.map { genome -> 
                                            def meta    = [:]
                                            meta.id     = genome.baseName
                                            
                                            [ meta, genome ]
                                        }

    // MODULE: SANITIZE_HEADERS
    SANITIZE_HEADERS ( ch_meta_genome )

    ch_sanitized_fasta                  = SANITIZE_HEADERS.out.fasta

    // MODULE: LTRHARVEST
    LTRHARVEST ( ch_sanitized_fasta )

    ch_ltrharvest_scn                   = LTRHARVEST.out.scn

    ch_versions                         = ch_versions.mix(LTRHARVEST.out.versions)

    // MODULE: LTRFINDER
    LTRFINDER  { ch_sanitized_fasta }

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
        ch_ltrretriever_inputs.map { meta, fasta, _ltr -> [ meta, fasta ] },
        ch_ltrretriever_inputs.map { _meta, _fasta, ltr -> ltr },
        [],
        [],
        []
    )

    ch_pass_list                        = LTRRETRIEVER_LTRRETRIEVER.out.pass_list
    ch_pass_list_gff                    = LTRRETRIEVER_LTRRETRIEVER.out.pass_list_gff
    ch_annotation_gff                   = LTRRETRIEVER_LTRRETRIEVER.out.annotation_gff
    ch_defalse                          = LTRRETRIEVER_LTRRETRIEVER.out.defalse
    ch_ltrlib                           = LTRRETRIEVER_LTRRETRIEVER.out.ltrlib
    ch_versions                         = ch_versions.mix(LTRRETRIEVER_LTRRETRIEVER.out.versions.first())

    // MODULE: TIRLEARNER
    TIRLEARNER (
        ch_sanitized_fasta,
        params.species
    )

    ch_tirlearner_fasta                 = TIRLEARNER.out.fasta
    ch_tirlearner_gff                   = TIRLEARNER.out.gff
    ch_versions                         = ch_versions.mix(TIRLEARNER.out.versions.first())

    // These can also run in parallel
    // MODULE: ANNOSINE
    ANNOSINE (
        ch_sanitized_fasta,
        3 // mode
    )

    // Currently it's a topic, so need to fix that
    ch_versions                         = ch_versions.mix(ANNOSINE.out.versions)
    cb_annosine_seed_sine               = ANNOSINE.out.fa

    // MODULE: REPEATMODELER_BUILDDATABASE
    ch_repeatmodeler_inputs             = ch_sanitized_fasta
                                        | map { meta, fasta ->
                                            def size = fasta.size()
                                            def size_threshold = 100_000 // bytes -> bp

                                            // TODO: Not the best way to set a size threshold
                                            // but it is simple
                                            // This is needed to avoid,
                                            // Error: Database genome is not large enough ( minimum 40000 bp ) to process with RepeatModeler.
                                            if ( size < size_threshold ) {
                                                log.warn "RepeatModeler is skipped for genome '${meta.id}' as it is smaller than ${size_threshold} bytes"
                                                return null
                                            }

                                            return [ meta, fasta ]
                                        }
                                        | filter { it }
    
    REPEATMODELER_BUILDDATABASE ( ch_repeatmodeler_inputs )

    ch_repeatmodeler_db                 = REPEATMODELER_BUILDDATABASE.out.db
    ch_versions                         = ch_versions.mix(REPEATMODELER_BUILDDATABASE.out.versions.first())

    // MODULE: REPEATMODELER_REPEATMODELER
    REPEATMODELER_REPEATMODELER ( ch_repeatmodeler_db )

    ch_repeatmodeler_fasta              = REPEATMODELER_REPEATMODELER.out.fasta
    ch_versions                         = ch_versions.mix(REPEATMODELER_REPEATMODELER.out.versions.first())

    // MODULE: FASTA_HELITRONSCANNER_SCAN_DRAW
    FASTA_HELITRONSCANNER_SCAN_DRAW ( ch_sanitized_fasta )

    ch_helitronscanner_draw             = FASTA_HELITRONSCANNER_SCAN_DRAW.out.helitronscanner_draw
    ch_helitronscanner_draw_rc          = FASTA_HELITRONSCANNER_SCAN_DRAW.out.helitronscanner_draw_rc
    ch_versions                         = ch_versions.mix(FASTA_HELITRONSCANNER_SCAN_DRAW.out.versions)

    // MODULE: FORMAT_HELITRONSCANNER_OUT
    ch_format_helitronscanner_inputs    = ch_sanitized_fasta
                                        | join(ch_helitronscanner_draw)
                                        | join(ch_helitronscanner_draw_rc)
                                        | multiMap { meta, fasta, draw, draw_rc ->
                                            genome: [ meta, fasta ]
                                            hel_fa: draw
                                            rc_hel_fa: draw_rc
                                        }
    
    FORMAT_HELITRONSCANNER_OUT (
        ch_format_helitronscanner_inputs.genome,
        ch_format_helitronscanner_inputs.hel_fa,
        ch_format_helitronscanner_inputs.rc_hel_fa,
    )

    ch_versions                         = ch_versions.mix(FORMAT_HELITRONSCANNER_OUT.out.versions.first())

    // MODULE: LTR_RETRIEVER_POSTPROCESS
    ch_ltr_retriever_postprocess_inputs = ch_sanitized_fasta
                                        | join(ch_pass_list)
                                        | join(ch_pass_list_gff)
                                        | join(ch_annotation_gff)
                                        | join(ch_defalse)
                                        | join(ch_ltrlib)
                                        | multiMap { meta, fasta, pass, p_gff, a_gff, defalse, ltr ->
                                            genome: [ meta, fasta ]
                                            pass    : pass
                                            p_gff   : p_gff
                                            a_gff   : a_gff
                                            defalse : defalse
                                            ltr     : ltr
                                        }
    
    LTR_RETRIEVER_POSTPROCESS (
        ch_ltr_retriever_postprocess_inputs.genome,
        ch_ltr_retriever_postprocess_inputs.pass,
        ch_ltr_retriever_postprocess_inputs.p_gff,
        ch_ltr_retriever_postprocess_inputs.defalse,
        ch_ltr_retriever_postprocess_inputs.ltr,
    )

    ch_versions                         = ch_versions.mix(LTR_RETRIEVER_POSTPROCESS.out.versions.first())

    // MODULE: TIR_LEARNER_POSTPROCESS
    ch_tir_learner_postprocess_inputs = ch_sanitized_fasta
                                        | join(ch_tirlearner_fasta)
                                        | join(ch_tirlearner_gff)
                                        | multiMap { meta, fasta, tir_fa, tir_gff ->
                                            genome  : [ meta, fasta ]
                                            tir_fa  : tir_fa
                                            tir_gff : tir_gff
                                        }
    TIR_LEARNER_POSTPROCESS (
        ch_tir_learner_postprocess_inputs.genome,
        ch_tir_learner_postprocess_inputs.tir_fa,
        ch_tir_learner_postprocess_inputs.tir_gff
    )

    ch_versions                         = ch_versions.mix(TIR_LEARNER_POSTPROCESS.out.versions.first())

    // Function: Save versions
    ch_versions                         = ch_versions
                                        | unique
                                        | map { yml ->
                                            if ( yml ) { yml }
                                        }
    
    ch_versions_yml                     = softwareVersionsToYAML(ch_versions)
                                        | collectFile(
                                            storeDir: "${params.outdir}/pipeline_info",
                                            name: 'software_versions.yml',
                                            sort: true,
                                            newLine: true,
                                            cache: false
                                        )

    emit:
    versions_yml                        = ch_versions_yml   // [ software_versions.yml ]

}
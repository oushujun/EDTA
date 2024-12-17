include { GUNZIP                                    } from '../modules/gallvp/gunzip/main'
include { CUSTOM_SHORTENFASTAIDS                    } from '../modules/gallvp/custom/shortenfastaids/main'
include { CUSTOM_RESTOREGFFIDS                      } from '../modules/gallvp/custom/restoregffids/main'

include { LTRHARVEST                                } from '../modules/nf-core/ltrharvest/main'
include { LTRFINDER                                 } from '../modules/nf-core/ltrfinder/main'
include { CAT_CAT                                   } from '../modules/nf-core/cat/cat/main'
include { LTRRETRIEVER_LTRRETRIEVER                 } from '../modules/nf-core/ltrretriever/ltrretriever/main'
include { LTR_RETRIEVER_POSTPROCESS                 } from '../modules/local/ltr_retriever_postprocess/main'

include { ANNOSINE                                  } from '../modules/gallvp/annosine/main'
include { ANNOSINE_POSTPROCESS                      } from '../modules/local/annosine_postprocess/main'

include { REPEATMODELER_BUILDDATABASE               } from '../modules/nf-core/repeatmodeler/builddatabase/main'
include { REPEATMODELER_REPEATMODELER               } from '../modules/nf-core/repeatmodeler/repeatmodeler/main'
include { REPEATMODELER_POSTPROCESS                 } from '../modules/local/repeatmodeler_postprocess/main'

include { TIRLEARNER                                } from '../modules/gallvp/tirlearner/main'
include { TIR_LEARNER_POSTPROCESS                   } from '../modules/local/tir_learner_postprocess/main'

include { FASTA_HELITRONSCANNER_SCAN_DRAW           } from '../subworkflows/gallvp/fasta_helitronscanner_scan_draw/main'
include { FORMAT_HELITRONSCANNER_OUT                } from '../modules/local/format_helitronscanner_out/main'
include { FORMAT_HELITRONSCANNER_OUT as FORMAT_HELITRONSCANNER_OUT_EXT  } from '../modules/local/format_helitronscanner_out/main'
include { HELITRONSCANNER_POSTPROCESS               } from '../modules/local/helitronscanner_postprocess/main'

include { COMBINE_INTACT_TES                        } from '../modules/local/combine_intact_tes/main'
include { PROCESS_K                                 } from '../modules/local/process_k/main'
include { FINAL_FILTER                              } from '../modules/local/final_filter/main'

include { softwareVersionsToYAML                    } from '../modules/local/utils/main'
include { idFromFileName                            } from '../modules/local/utils/main'

workflow EDTA {

    main:

    // Versions channel
    ch_versions                                     = Channel.empty()

    
    ch_genome_branch                                = Channel.fromPath(params.genome)
                                                    | map { genome -> 
                                                        def meta    = [:]
                                                        meta.id     = idFromFileName ( genome.baseName )
                                                        
                                                        [ meta, genome ]
                                                    }
                                                    | branch { _meta, archive ->
                                                        gz: "$archive".endsWith('.gz')
                                                        rest: ! "$archive".endsWith('.gz')
                                                    }

    // MODULE: GUNZIP                   
    GUNZIP ( ch_genome_branch.gz )

    ch_genome               = GUNZIP.out.gunzip.mix(ch_genome_branch.rest)
    ch_versions             = ch_versions.mix(GUNZIP.out.versions.first())

    // MODULE: CUSTOM_SHORTENFASTAIDS
    CUSTOM_SHORTENFASTAIDS ( ch_genome )

    ch_short_ids_tsv                                = CUSTOM_SHORTENFASTAIDS.out.short_ids_tsv
    ch_versions                                     = ch_versions.mix(CUSTOM_SHORTENFASTAIDS.out.versions.first())

    
    ch_shortenfastaids_branch                       = ch_short_ids_tsv
                                                    | branch { _meta, tsv ->
                                                        change: ! tsv.text.contains('IDs have acceptable length and character')
                                                        nochange: tsv.text.contains('IDs have acceptable length and character')
                                                    }

    ch_sanitized_fasta                              = ch_shortenfastaids_branch.nochange
                                                    | join(
                                                        ch_genome
                                                    )
                                                    | map { meta, _tsv, fasta -> [ meta + [ changed_ids: false ], fasta ] }
                                                    | mix(
                                                        ch_shortenfastaids_branch.change
                                                        | join(
                                                            CUSTOM_SHORTENFASTAIDS.out.short_ids_fasta
                                                        )
                                                        | map { meta, _tsv, fasta -> [ meta + [ changed_ids: true ], fasta ] }
                                                    )

    // MODULE: LTRHARVEST
    LTRHARVEST ( ch_sanitized_fasta )

    ch_ltrharvest_scn                               = LTRHARVEST.out.scn

    ch_versions                                     = ch_versions.mix(LTRHARVEST.out.versions.first())

    // MODULE: LTRFINDER
    LTRFINDER  { ch_sanitized_fasta }

    ch_ltrfinder_scn                                = LTRFINDER.out.scn

    ch_versions                                     = ch_versions.mix(LTRFINDER.out.versions.first())

    // MODULE: CAT_CAT
    ch_cat_cat_inputs                               = ch_ltrharvest_scn
                                                    | join(ch_ltrfinder_scn)
                                                    | map { meta, harvested, found -> [ meta, [ harvested, found ] ] }
    CAT_CAT ( ch_cat_cat_inputs )

    ch_ltr_candidates                               = CAT_CAT.out.file_out
    ch_versions                                     = ch_versions.mix(CAT_CAT.out.versions.first())

    // MODULE: LTRRETRIEVER_LTRRETRIEVER
    ch_ltrretriever_inputs                          = ch_sanitized_fasta.join(ch_ltr_candidates)

    LTRRETRIEVER_LTRRETRIEVER (
        ch_ltrretriever_inputs.map { meta, fasta, _ltr -> [ meta, fasta ] },
        ch_ltrretriever_inputs.map { _meta, _fasta, ltr -> ltr },
        [],
        [],
        []
    )

    ch_versions                                     = ch_versions.mix(LTRRETRIEVER_LTRRETRIEVER.out.versions.first())

    // MODULE: LTR_RETRIEVER_POSTPROCESS
    ch_ltr_retriever_postprocess_inputs             = ch_sanitized_fasta
                                                    | join ( LTRRETRIEVER_LTRRETRIEVER.out.pass_list        )
                                                    | join ( LTRRETRIEVER_LTRRETRIEVER.out.pass_list_gff    )
                                                    | join ( LTRRETRIEVER_LTRRETRIEVER.out.annotation_gff   )
                                                    | join ( LTRRETRIEVER_LTRRETRIEVER.out.defalse          )
                                                    | join ( LTRRETRIEVER_LTRRETRIEVER.out.ltrlib           )
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

    ch_versions                                     = ch_versions.mix(LTR_RETRIEVER_POSTPROCESS.out.versions.first())

    // MODULE: ANNOSINE
    ANNOSINE (
        ch_sanitized_fasta,
        3 // mode
    )

    ch_versions                                     = ch_versions.mix(ANNOSINE.out.versions.first())

    // MODULE: ANNOSINE_POSTPROCESS
    ch_annosine_postprocess_inputs                  = ch_sanitized_fasta
                                                    | join ( ANNOSINE.out.fa )
                                                    | multiMap { meta, fasta, anno_fa ->
                                                        genome  : [ meta, fasta ]
                                                        anno_fa : anno_fa
                                                    }
    ANNOSINE_POSTPROCESS (
        ch_annosine_postprocess_inputs.genome,
        ch_annosine_postprocess_inputs.anno_fa
    )

    ch_versions                                     = ch_versions.mix(ANNOSINE_POSTPROCESS.out.versions.first())

    // MODULE: REPEATMODELER_BUILDDATABASE
    ch_repeatmodeler_inputs                         = ch_sanitized_fasta
                                                    | map { meta, fasta ->
                                                        def size = fasta.size()
                                                        def size_threshold = 100_000 // bytes -> bp

                                                        if ( size < size_threshold ) {
                                                            log.warn "RepeatModeler is skipped for genome '${meta.id}' as it is smaller than ${size_threshold} bytes"
                                                            return null
                                                        }

                                                        return [ meta, fasta ]
                                                    }
                                                    | filter { it }
    
    REPEATMODELER_BUILDDATABASE ( ch_repeatmodeler_inputs )

    ch_repeatmodeler_db                             = REPEATMODELER_BUILDDATABASE.out.db
    ch_versions                                     = ch_versions.mix(REPEATMODELER_BUILDDATABASE.out.versions.first())

    // MODULE: REPEATMODELER_REPEATMODELER
    REPEATMODELER_REPEATMODELER ( ch_repeatmodeler_db )

    ch_versions                                     = ch_versions.mix(REPEATMODELER_REPEATMODELER.out.versions.first())

    // MODULE: REPEATMODELER_POSTPROCESS
    ch_repeatmodeler_postprocess_inputs             = ch_sanitized_fasta
                                                    | join ( REPEATMODELER_REPEATMODELER.out.fasta )
                                                    | multiMap { meta, fasta, rm_fa ->
                                                        genome  : [ meta, fasta ]
                                                        rm_fa   : rm_fa
                                                    }
    REPEATMODELER_POSTPROCESS (
        ch_repeatmodeler_postprocess_inputs.genome,
        ch_repeatmodeler_postprocess_inputs.rm_fa,
    )

    ch_versions                                     = ch_versions.mix(REPEATMODELER_POSTPROCESS.out.versions.first())

    // MODULE: TIRLEARNER
    TIRLEARNER (
        ch_sanitized_fasta,
        params.species
    )

    ch_versions                                     = ch_versions.mix(TIRLEARNER.out.versions.first())

    // MODULE: TIR_LEARNER_POSTPROCESS
    ch_tir_learner_postprocess_inputs               = ch_sanitized_fasta
                                                    | join ( TIRLEARNER.out.fasta   )
                                                    | join ( TIRLEARNER.out.gff     )
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

    ch_versions                                     = ch_versions.mix(TIR_LEARNER_POSTPROCESS.out.versions.first())

    // MODULE: FASTA_HELITRONSCANNER_SCAN_DRAW
    FASTA_HELITRONSCANNER_SCAN_DRAW ( ch_sanitized_fasta )

    ch_helitronscanner_draw                         = FASTA_HELITRONSCANNER_SCAN_DRAW.out.helitronscanner_draw
    ch_helitronscanner_draw_rc                      = FASTA_HELITRONSCANNER_SCAN_DRAW.out.helitronscanner_draw_rc
    ch_versions                                     = ch_versions.mix(FASTA_HELITRONSCANNER_SCAN_DRAW.out.versions)

    // MODULE: FORMAT_HELITRONSCANNER_OUT
    ch_format_helitronscanner_inputs                = ch_sanitized_fasta
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

    ch_helitronscanner_out_fa                       = FORMAT_HELITRONSCANNER_OUT.out.filtered_fa
                                                    | filter { _meta, fasta -> fasta.countFasta() > 0 }
    ch_versions                                     = ch_versions.mix(FORMAT_HELITRONSCANNER_OUT.out.versions.first())
    
    
    // MODULE: FORMAT_HELITRONSCANNER_OUT as FORMAT_HELITRONSCANNER_OUT_EXT
    FORMAT_HELITRONSCANNER_OUT_EXT (
        ch_format_helitronscanner_inputs.genome,
        ch_format_helitronscanner_inputs.hel_fa,
        ch_format_helitronscanner_inputs.rc_hel_fa,
    )

    ch_helitronscanner_out_ext_fa                   = FORMAT_HELITRONSCANNER_OUT_EXT.out.filtered_ext_fa
                                                    | filter { _meta, fasta -> fasta.countFasta() > 0 }
    ch_versions                                     = ch_versions.mix(FORMAT_HELITRONSCANNER_OUT_EXT.out.versions.first())

    // MODULE: HELITRONSCANNER_POSTPROCESS
    ch_helitronscanner_post_inputs                  = ch_sanitized_fasta
                                                    | join(ch_helitronscanner_out_fa)
                                                    | join(ch_helitronscanner_out_ext_fa)
                                                    | multiMap { meta, fasta, hs_fa, hs_ext_fa ->
                                                        genome      : [ meta, fasta ]
                                                        hs_fa       : hs_fa
                                                        hs_ext_fa   : hs_ext_fa
                                                    }
    HELITRONSCANNER_POSTPROCESS (
        ch_helitronscanner_post_inputs.genome,
        ch_helitronscanner_post_inputs.hs_fa,
        ch_helitronscanner_post_inputs.hs_ext_fa,
    )

    ch_versions                                     = ch_versions.mix(HELITRONSCANNER_POSTPROCESS.out.versions.first())

    // MODULE: COMBINE_INTACT_TES
    ch_combine_intact_tes_inputs                    = ch_sanitized_fasta
                                                    | join ( LTR_RETRIEVER_POSTPROCESS.out.intact_raw_fa    )
                                                    | join ( LTR_RETRIEVER_POSTPROCESS.out.intact_raw_gff3  )
                                                    | join ( TIR_LEARNER_POSTPROCESS.out.intact_raw_fa      )
                                                    | join ( TIR_LEARNER_POSTPROCESS.out.intact_raw_bed     )
                                                    | join ( HELITRONSCANNER_POSTPROCESS.out.raw_fa         )
                                                    | join ( HELITRONSCANNER_POSTPROCESS.out.raw_bed        )
                                                    | multiMap { meta, genome, ltr_fa, ltr_gff, tir_fa, tir_bed, helitron_fa, helitron_bed ->
                                                        genome      : [ meta, genome ]
                                                        ltr_fa      : ltr_fa
                                                        ltr_gff     : ltr_gff
                                                        tir_fa      : tir_fa
                                                        tir_bed     : tir_bed
                                                        helitron_fa : helitron_fa
                                                        helitron_bed: helitron_bed
                                                    }

    COMBINE_INTACT_TES (
        ch_combine_intact_tes_inputs.genome,
        ch_combine_intact_tes_inputs.ltr_fa,
        ch_combine_intact_tes_inputs.ltr_gff,
        ch_combine_intact_tes_inputs.tir_fa,
        ch_combine_intact_tes_inputs.tir_bed,
        ch_combine_intact_tes_inputs.helitron_fa,
        ch_combine_intact_tes_inputs.helitron_bed,
    )

    ch_versions                                     = ch_versions.mix(COMBINE_INTACT_TES.out.versions.first())

    // MODULE: PROCESS_K
    ch_process_k_inputs                             = ch_sanitized_fasta
                                                    | join ( LTR_RETRIEVER_POSTPROCESS.out.raw_fa                           )
                                                    | join ( LTR_RETRIEVER_POSTPROCESS.out.intact_raw_fa                    )
                                                    | join ( ANNOSINE_POSTPROCESS.out.sine_fa           , remainder: true   )
                                                    | join ( REPEATMODELER_POSTPROCESS.out.line_raw     , remainder: true   )
                                                    | join ( TIR_LEARNER_POSTPROCESS.out.intact_raw_fa                      )
                                                    | join ( HELITRONSCANNER_POSTPROCESS.out.raw_fa                         )
                                                    | multiMap { meta, genome, ltr, ltrint, sine, line, tir, helitron ->
                                                        genome  : [ meta, genome ]
                                                        ltr     : ltr       ?: []
                                                        ltrint  : ltrint    ?: []
                                                        sine    : sine      ?: []
                                                        line    : line      ?: []
                                                        tir     : tir       ?: []
                                                        helitron: helitron  ?: []
                                                    }

    PROCESS_K (
        ch_process_k_inputs.genome,
        ch_process_k_inputs.ltr,
        ch_process_k_inputs.ltrint,
        ch_process_k_inputs.sine,
        ch_process_k_inputs.line,
        ch_process_k_inputs.tir,
        ch_process_k_inputs.helitron,
    )

    ch_versions                                     = ch_versions.mix(PROCESS_K.out.versions.first())

    // Warn: End of the pipeline if no inputs are available for ProcessK
    ch_sanitized_fasta
    | join ( LTR_RETRIEVER_POSTPROCESS.out.raw_fa           , remainder: true )
    | join ( LTR_RETRIEVER_POSTPROCESS.out.intact_raw_fa    , remainder: true )
    | join ( TIR_LEARNER_POSTPROCESS.out.intact_raw_fa      , remainder: true )
    | join ( HELITRONSCANNER_POSTPROCESS.out.raw_fa         , remainder: true )
    | map { meta, _genome, ltr, ltrint, tir, helitron ->
        if ( !ltr || !ltrint || !tir || !helitron ) {
            log.warn "One or more TE classes needed to complete EDTA were not found in genome '$meta.id'" +
            ". Multiple processing steps are being skipped."
        }
    }

    // MODULE: FINAL_FILTER
    ch_final_filter_inputs                          = ch_sanitized_fasta
                                                    | join(PROCESS_K.out.stg1_fa)
                                                    | join(PROCESS_K.out.intact_fa)
                                                    | join(COMBINE_INTACT_TES.out.intact_raw_gff)
                                                    | multiMap { meta, fasta, stg1_fa, intact_fa, intact_gff ->
                                                        genome      : [ meta, fasta ]
                                                        stg1_fa     : stg1_fa
                                                        intact_fa   : intact_fa
                                                        intact_gff  : intact_gff
                                                    }
    FINAL_FILTER (
        ch_final_filter_inputs.genome,
        ch_final_filter_inputs.stg1_fa,
        ch_final_filter_inputs.intact_fa,
        ch_final_filter_inputs.intact_gff,
    )

    ch_intact_gff                                   = FINAL_FILTER.out.intact_gff
    ch_versions                                     = ch_versions.mix(FINAL_FILTER.out.versions.first())

    // MODULE: CUSTOM_RESTOREGFFIDS
    ch_gff_tsv_branch               = ch_intact_gff.join(ch_short_ids_tsv)
                                    | branch { meta, _gff, _tsv ->
                                        change: meta.changed_ids
                                        nochange: ! meta.changed_ids
                                    }

    CUSTOM_RESTOREGFFIDS (
        ch_gff_tsv_branch.change.map { meta, gff, _tsv -> [ meta, gff ] },
        ch_gff_tsv_branch.change.map { _meta, _gff, tsv -> tsv }
    )

    ch_versions                     = ch_versions.mix(CUSTOM_RESTOREGFFIDS.out.versions.first())

    // Function: Save versions
    ch_versions                                     = ch_versions
                                                    | unique
                                                    | map { yml ->
                                                        if ( yml ) { yml }
                                                    }
    
    ch_versions_yml                                 = softwareVersionsToYAML(ch_versions)
                                                    | collectFile(
                                                        storeDir: "${params.outdir}/pipeline_info",
                                                        name: 'software_versions.yml',
                                                        sort: true,
                                                        newLine: true,
                                                        cache: false
                                                    )

    emit:
    versions_yml                                    = ch_versions_yml   // [ software_versions.yml ]

}
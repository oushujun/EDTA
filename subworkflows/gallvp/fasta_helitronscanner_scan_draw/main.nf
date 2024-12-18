include { HELITRONSCANNER_SCAN as HELITRONSCANNER_SCAN_HEAD     } from '../../../modules/gallvp/helitronscanner/scan/main.nf'
include { HELITRONSCANNER_SCAN as HELITRONSCANNER_SCAN_TAIL     } from '../../../modules/gallvp/helitronscanner/scan/main.nf'
include { HELITRONSCANNER_DRAW                                  } from '../../../modules/gallvp/helitronscanner/draw/main.nf'

include { HELITRONSCANNER_SCAN as HELITRONSCANNER_SCAN_HEAD_RC  } from '../../../modules/gallvp/helitronscanner/scan/main.nf'
include { HELITRONSCANNER_SCAN as HELITRONSCANNER_SCAN_TAIL_RC  } from '../../../modules/gallvp/helitronscanner/scan/main.nf'
include { HELITRONSCANNER_DRAW as HELITRONSCANNER_DRAW_RC       } from '../../../modules/gallvp/helitronscanner/draw/main.nf'

workflow FASTA_HELITRONSCANNER_SCAN_DRAW {

    take:
    ch_fasta                        // channel: [ val(meta), fasta ]

    main:

    ch_versions = Channel.empty()

    // MODULE: HELITRONSCANNER_SCAN  as HELITRONSCANNER_SCAN_HEAD
    HELITRONSCANNER_SCAN_HEAD (
        ch_fasta,
        'head',     // command
        [],         // lcv_filepath
        0           // buffer_size
    )

    ch_helitronscanner_scan_head    = HELITRONSCANNER_SCAN_HEAD.out.scan
    ch_versions                     = ch_versions.mix(HELITRONSCANNER_SCAN_HEAD.out.versions)

    // MODULE: HELITRONSCANNER_SCAN  as HELITRONSCANNER_SCAN_TAIL
    HELITRONSCANNER_SCAN_TAIL (
        ch_fasta,
        'tail',     // command
        [],         // lcv_filepath
        0           // buffer_size
    )

    ch_helitronscanner_scan_tail    = HELITRONSCANNER_SCAN_TAIL.out.scan
    ch_versions                     = ch_versions.mix(HELITRONSCANNER_SCAN_TAIL.out.versions)

    // MODULE: HELITRONSCANNER_DRAW
    HELITRONSCANNER_DRAW (
        ch_fasta,
        ch_helitronscanner_scan_head,
        ch_helitronscanner_scan_tail
    )

    ch_helitronscanner_draw         = HELITRONSCANNER_DRAW.out.draw
    ch_versions                     = ch_versions.mix(HELITRONSCANNER_DRAW.out.versions)

    // MODULE: HELITRONSCANNER_SCAN  as HELITRONSCANNER_SCAN_HEAD_RC
    HELITRONSCANNER_SCAN_HEAD_RC (
        ch_fasta,
        'head',     // command
        [],         // lcv_filepath
        0           // buffer_size
    )

    ch_helitronscanner_scan_head_rc = HELITRONSCANNER_SCAN_HEAD_RC.out.scan
    ch_versions                     = ch_versions.mix(HELITRONSCANNER_SCAN_HEAD_RC.out.versions)

    // MODULE: HELITRONSCANNER_SCAN  as HELITRONSCANNER_SCAN_TAIL_RC
    HELITRONSCANNER_SCAN_TAIL_RC (
        ch_fasta,
        'tail',     // command
        [],         // lcv_filepath
        0           // buffer_size
    )

    ch_helitronscanner_scan_tail_rc = HELITRONSCANNER_SCAN_TAIL_RC.out.scan
    ch_versions                     = ch_versions.mix(HELITRONSCANNER_SCAN_TAIL_RC.out.versions)

    // MODULE: HELITRONSCANNER_DRAW as HELITRONSCANNER_DRAW_RC
    HELITRONSCANNER_DRAW_RC (
        ch_fasta,
        ch_helitronscanner_scan_head_rc,
        ch_helitronscanner_scan_tail_rc
    )

    ch_helitronscanner_draw_rc      = HELITRONSCANNER_DRAW_RC.out.draw
    ch_versions                     = ch_versions.mix(HELITRONSCANNER_DRAW_RC.out.versions)

    emit:
    helitronscanner_draw            = ch_helitronscanner_draw       // channel: [ val(meta), draw ]
    helitronscanner_draw_rc         = ch_helitronscanner_draw_rc    // channel: [ val(meta), rc.draw ]
    versions                        = ch_versions                   // channel: [ versions.yml ]
}

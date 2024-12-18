include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

workflow PIPELINE_INITIALISATION {

    take:
    version                 // boolean: Display version and exit
    nextflow_cli_args       //   array: List of positional nextflow CLI args
    outdir                  //  string: The output directory where the results will be saved
    genome                  //  string: Path to a genome fasta or a text file with a list of fasta files

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        true, // validate params
        null
    )

    //
    // Create input channels
    //
    ch_genome                   = genome.find(/^\S+\.f(a|asta|as|sa|na)(\.gz)?$/)
                                ? (
                                    // Fasta
                                    Channel.fromPath(genome)
                                    | map { fasta -> 
                                        def fa_file = file(fasta, checkIfExists: true)
                                        
                                        def meta    = [:]
                                        meta.id     = idFromFileName ( fa_file.baseName )
                                        
                                        [ meta, fa_file ]
                                    }
                                )
                                : (
                                    // Text file
                                    Channel.fromPath(genome)
                                    | map { text ->
                                        def file_list = file(text, checkIfExists: true).readLines()
                                            
                                        if ( file_list == [] ) {
                                            error "File list provided by --genome is empty"
                                        }
                                        
                                        file_list
                                            .collect { fasta ->
                                                def fa_file = file(fasta, checkIfExists: true)
                                                
                                                def meta    = [:]
                                                meta.id     = idFromFileName ( fa_file.baseName )
                                        
                                                [ meta, fa_file ]
                                            }
                                    }
                                    | flatten
                                    | buffer ( size: 2 )
                                )
                                | groupTuple
                                | view
                                | map { meta, fastas -> validateFastaMetadata ( meta, fastas ) }
    
    emit:
    genome                      = ch_genome
}

//
// Additional validation
//
def idFromFileName(fileName) {

    def trial = ( fileName
        ).replaceFirst(
            /\.f(ast)?q$/, ''
        ).replaceFirst(
            /\.f(asta|sa|a|as|aa|na)?$/, ''
        ).replaceFirst(
            /\.gff(3)?$/, ''
        ).replaceFirst(
            /\.gz$/, ''
        )

    if ( trial == fileName ) { return fileName }

    return idFromFileName ( trial )
}

def validateFastaMetadata(meta, fastas) {

    if ( fastas.size() > 1 ) {
        error "Multiple Fasta files have the same ID '${meta.id}'. Make sure all the files are unique and change file names to avoid this conflict!"
    }

    [ meta, fastas.first() ]
}

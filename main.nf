!/usr/bin/env nextflow

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

// TODO: Check parameters, also help
// - We should at a later stage use [nf-validation](https://github.com/nextflow-io/nf-validation) and
// nf-core pipeline [template](https://nf-co.re/docs/contributing/adding_pipelines#create-a-pipeline-from-the-template)
// for parameter validation.
//
// - I have implemented this for [plant-food-research-open/assemblyqc](https://github.com/Plant-Food-Research-Open/assemblyqc)
//
// - I think we should also aspire to support docker, conda and singularity engines
//
// - We should follow the nf-core module pattern, e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/ltrretriever/ltrretriever/main.nf
/*

if ($maxdiv < 0 or $maxdiv > 100){die "The expected value for the div parameter is 0 - 100!\n"}
if ($overwrite != 0 and $overwrite != 1){ die "The expected value for the overwrite parameter is 0 or 1!\n"}
if ($sensitive != 0 and $sensitive != 1){ die "The expected value for the sensitive parameter is 0 or 1!\n"}
if ($anno != 0 and $anno != 1){ die "The expected value for the anno parameter is 0 or 1!\n"}
if ($evaluate != 0 and $evaluate != 1){ die "The expected value for the evaluate parameter is 0 or 1!\n"}
if ($force != 0 and $force != 1){ die "The expected value for the force parameter is 0 or 1!\n"}
if ($miu !~ /[0-9\.e\-]+/){ die "The expected value for the u parameter is float value without units!\n"}
if ($debug != 0 and $debug != 1){ die "The expected value for the debug parameter is 0 or 1!\n"}
if ($threads !~ /^[0-9]+$/){ die "The expected value for the threads parameter is an integer!\n"}

*/

// TODO: Check inputed repeat libraries, CDS, etc...
// TODO: Check exclude file

// Rename FASTA headers (just makes everything easier later)
// TODO: Put fffx on bioconda or somewhere so it just runs, otherwise tiny container
// NEW: Filter out small contigs (1000bp or less, remove)
// TODO: Make sanitize command do filtering, bbmap is memory hungry!
// fffx is compiled for Linux x86_64? I'll be developing on macos amd64.
// We'll have to put the the software in conda/docker containers
process sanitize {
    tag "${x.baseName}"
    input:
        path x
    output:
        tuple val("${x.baseName}"), path("${x.baseName}_sanitized.fasta"), path("${x.baseName}_sanitized.translation_table.tsv")
    publishDir 'sanitized_genomes'
    time "10m"
    memory 3.GB
    cpus 1

"""
fffx length-filter ${x} filtered.fa 1000
fffx sanitize filtered.fa ${x.baseName}_sanitized
"""
}

// It's actually single threaded
process helitron_scanner {
    tag "${genome.baseName}"
    input:
        path(genome)
    output:
        path("${genome.baseName}.Helitron.intact.raw.gff3")
    conda 'bioconda::tesorter bioconda::mdust bioconda::trf'
    cpus 1 
    time '18h'
    memory 32.GB
    publishDir 'out_helitron_scanner'

"""
sh ${projectDir}/util/run_helitron_scanner.sh \
    ${genome} \
    ${task.cpus}

perl ${projectDir}/util/format_helitronscanner_out.pl \
    -genome $genome \
    -sitefilter 1 \
    -minscore 12 \
    -keepshorter 1 \
    -extlen 30 \
    -extout 1

perl ${projectDir}/util/flanking_filter.pl \
    -genome ${genome} \
    -query ${genome}.HelitronScanner.filtered.ext.fa \
    -miniden 90 \
    -mincov 0.9 \
    -maxct 5 \
    -t ${task.cpus}

# remove simple repeats and candidates with simple repeats at terminals
perl ${projectDir}/util/output_by_list.pl 1 \
    ${genome}.HelitronScanner.filtered.fa \
    1 \
    ${genome}.HelitronScanner.filtered.ext.fa \
    -FA > ${genome}.HelitronScanner.filtered.pass.fa

mdust ${genome}.HelitronScanner.filtered.pass.fa > ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted
perl ${projectDir}/util/cleanup_tandem.pl \
    -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 100 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted | perl ${projectDir}/util/helitron_renamer.pl > ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln
# Too complicated for nextflow
# perl -nle 's/^(>.*)\\s+(.*)\\\$/\\\$1#DNA\\/Helitron\\t\\\$2/; print \\\$_' > ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln

# annotate and remove non-Helitron candidates
TEsorter ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln \
    --disable-pass2 \
    -p ${task.cpus}

perl ${projectDir}/util/cleanup_misclas.pl \
    ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln.rexdb.cls.tsv
mv ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln ${genome.baseName}.Helitron.intact.raw.fa
cp ${genome}.HelitronScanner.filtered.fa.pass.fa.dusted.cln.cln.list ${genome.baseName}.Helitron.intact.raw.fa.anno.list

# get intact Helitrons and gff3
perl ${projectDir}/util/make_bed_with_intact.pl \
    ${genome.baseName}.Helitron.intact.raw.fa > ${genome.baseName}.Helitron.intact.raw.bed
perl ${projectDir}/util/bed2gff.pl ${genome.baseName}.Helitron.intact.raw.bed HEL > ${genome.baseName}.Helitron.intact.raw.gff3
"""
}

// Includes
include { LTR_HARVEST                   } from './subworkflows/local/ltr_harvest'
include { LTR_FINDER                    } from './subworkflows/local/ltr_finder'
include { LTR_RETRIEVER                 } from './subworkflows/local/ltr_retriever'
include { ANNOSINE                      } from './subworkflows/local/annosine'
include { REPEATMODELER                 } from './subworkflows/local/repeatmodeler'
include { TIRLEARNER                    } from './subworkflows/local/tir_learner'

workflow WORKFLOW_A {
    // - All nf-core pipelines/modules the [ [`meta`], data] pattern. I think we should follow that
    // so that our local and nf-core modules are interoperable
    //
    // - I am also adding a bit of personal code style, please rever it if you don't like it

    ch_genome                           = Channel.fromPath(params.genomes) // The channel emits a single genome at a time. That's why ch_genome
    
    // MODULE: sanitize - All modules and workflow names should follow CAPITAL_SNAKE
    sanitize ( ch_genome )
 

    // Handy for data control
    ch_genome_all                       = sanitize.out
                                        | map { id, fasta, table ->
                                            [
                                                id,
                                                [ "name": id, "assembly": fasta, "translation_table": table ],
                                                fasta
                                            ]
                                        }

    // ch_genomes_input                    = sanitize.out
    //                                     | map({ it ->
    //     tuple(["name": it[0], "assembly": it[1], "translation_table": it[2]], it[1])
    // }) // Smame as ch_genomes_all?

    ch_fasta                                = sanitize.out
                                            | map { id, fasta, table ->
                                                [ id, fasta ]
                                            }
    
    // Input for (nearly all) functions are going to be the map "data"
    // Output will be a tuple (data.name, output_file)
    // So we can merge them all later using data.name

    // The map has the folliwing features:
    // - name: The name of the genome
    // - assembly: The sanitized genome
    // - translation_table: The translation table for the genome

    // These two feed into ltr_retriever (but the first two run in parallel)
    // SUBWORKFLOW: LTR_HARVEST
    LTR_HARVEST ( ch_genome_all )

    //
    // AT THIS POINT
    // Only go to this portion for development


    // SUBWORKFLOW: LTR_FINDER
    LTR_FINDER ( ch_genome_all )
    
    ch_ltr_retriever_input                  = LTR_HARVEST.out
                                            | join(LTR_FINDER.out)
                                            | join(ch_fasta)
    
    // SUBWORKFLOW: LTR_RETRIEVER
    LTR_RETRIEVER ( ch_ltr_retriever_input )

    // These can also run in parallel
    // SUBWORKFLOW: ANNOSINE
    ANNOSINE (ch_genome_all )
    
    // SUBWORKFLOW: REPEATMODELER
    REPEATMODELER ( ch_genome_all )
    
    // SUBWORKFLOW: TIRLEARNER
    TIRLEARNER ( ch_genome_all )
    // helitron_scanner(sanitized_genomes)
}

// Functions
include { idFromFileName                    } from './modules/local/utils'

// MODULES
include { GUNZIP                            } from './modules/nf-core/gunzip/main'
include { CUSTOM_SHORTENFASTAIDS            } from './modules/pfr/custom/shortenfastaids/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from './modules/nf-core/custom/dumpsoftwareversions/main'

// SUBWORKFLOWS
include { FASTA_LTRRETRIEVER                } from './subworkflows/local/fasta_ltrretriever'
include { FASTA_REPEATMODELER               } from './subworkflows/local/fasta_repeatmodeler'
include { FASTA_ANNOSINE                    } from './subworkflows/local/fasta_annosine'

workflow WORKFLOW_B {
    // - Using nf-core styling and modules which support conda, docker and singularity
    // 
    // - To test it with docker run: ./main.nf -profile local,docker -resume -c conf/test.config
    // 
    // - To test it with conda run: ./main.nf -profile local,conda -resume -c conf/test.config
    //
    // - Tools installed by nf-core -v modules install <tool/subtool>
    //
    // - stub is also supported
    // Test: ./main.nf -profile local,docker -resume -stub -c conf/test.config

    // Versions
    ch_versions                             = Channel.empty()
    
    // Input channels
    ch_fasta_input                          = Channel.fromPath(params.genomes)
                                            | map { fasta ->
                                                [ [ id: idFromFileName(fasta.baseName) ], fasta ]
                                            }

    ch_fasta_branch                         = ch_fasta_input
                                            | branch { meta, fasta ->
                                                gz: "$fasta".endsWith(".gz")
                                                rest: ! "$fasta".endsWith(".gz")
                                            }

    // MODULE: GUNZIP
    GUNZIP ( ch_fasta_branch.gz )

    ch_fasta                                = GUNZIP.out.gunzip.mix(ch_fasta_branch.rest)
    ch_versions                             = ch_versions.mix(GUNZIP.out.versions.first())

    // MOUDLE: CUSTOM_SHORTENFASTAIDS
    CUSTOM_SHORTENFASTAIDS ( ch_fasta )

    ch_short_ids_fasta                      = ch_fasta
                                            | join(CUSTOM_SHORTENFASTAIDS.out.short_ids_fasta, by:0, remainder:true)
                                            | map { meta, fasta, short_ids_fasta ->
                                                if ( fasta ) { [ meta, short_ids_fasta ?: fasta ] }
                                            }

    ch_short_ids_tsv                        = CUSTOM_SHORTENFASTAIDS.out.short_ids_tsv
    ch_versions                             = ch_versions.mix(CUSTOM_SHORTENFASTAIDS.out.versions.first())

    // SUBWORKFLOW: FASTA_LTRRETRIEVER
    FASTA_LTRRETRIEVER ( ch_short_ids_fasta, ch_short_ids_tsv )

    ch_versions                             = ch_versions.mix(FASTA_LTRRETRIEVER.out.versions)

    // SUBWORKFLOW: FASTA_REPEATMODELER
    FASTA_REPEATMODELER ( ch_short_ids_fasta )

    ch_versions                             = ch_versions.mix(FASTA_REPEATMODELER.out.versions)

    // SUBWORKFLOW: FASTA_ANNOSINE
    FASTA_ANNOSINE ( ch_short_ids_fasta )

    ch_versions                             = ch_versions.mix(FASTA_ANNOSINE.out.versions)

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

workflow {
    // WORKFLOW_A ( ) // Uncomment to use WORKFLOW_A
    WORKFLOW_B ( )
}
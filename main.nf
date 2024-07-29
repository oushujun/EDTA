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

// TODO: Check inputed repeat libraries, CDS, etc...
// TODO: Check exclude file

// Rename FASTA headers (just makes everything easier later)
// TODO: Put fffx on bioconda or somewhere so it just runs, otherwise tiny container
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

workflow WORKFLOW_A {
    // - All nf-core pipelines/modules the [ [`meta`], data] pattern. I think we should follow that
    // so that our local and nf-core modules are interoperable
    //
    // - I am also adding a bit of personal code style, please rever it if you don't like it

    ch_genome                           = Channel.fromPath(params.genomes) // The channel emits a single genome at a time. That's why ch_genome
    
    // MODULE: sanitize - All modules and workflow names should follow CAPITAL_SNAKE
    sanitize ( ch_genome )
}
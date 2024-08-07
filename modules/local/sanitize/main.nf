// Rename FASTA headers (just makes everything easier later)

process SANITIZE_HEADERS {
    tag "$meta.id"
    label 'process_single'

    // Eventually port fffx (pronounced f3x) to bioconda
    // conda "${moduleDir}/environment.yml"
    // container "docker.io/gallvp/edta-components:v0.1"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04' }"
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.sanitized.fasta'), emit: fasta
    tuple val(meta), path('*.sanitized.translation_table.tsv'), emit: translation_table
    eval('fffx --version'), topic: versions 

    // todo use nf-core resource profiles?
    time "10m"
    memory 3.GB
    cpus 1

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    fffx length-filter ${fasta} filtered.fa 1000
    fffx sanitize filtered.fa ${fasta.baseName}.sanitized
    """

    stub:
    """
    touch ${fasta.baseName}.sanitized.fasta
    touch ${fasta.baseName}.sanitized.translation_table.tsv
    """
}
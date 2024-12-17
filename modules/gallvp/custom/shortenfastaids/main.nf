process CUSTOM_SHORTENFASTAIDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.short.ids.fasta")  , emit: short_ids_fasta , optional: true
    tuple val(meta), path("*.short.ids.tsv")    , emit: short_ids_tsv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'shorten_fasta_ids.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo \\
        'IDs have acceptable length and character. No change required.' \\
        > ${meta.id}.short.ids.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        biopython: \$(pip list | grep "biopython" | cut -d' ' -f3)
    END_VERSIONS
    """
}

process CUSTOM_RESTOREGFFIDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.2':
        'biocontainers/python:3.10.2' }"

    input:
    tuple val(meta), path(gff3)
    path(ids_tsv)

    output:
    tuple val(meta), path("*.restored.ids.gff3")    , emit: restored_ids_gff3
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'restore_gff_ids.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.restored.ids.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}

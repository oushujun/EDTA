process REPEATMODELER_BUILDDATABASE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.5--pl5321hdfd78af_0':
        'biocontainers/repeatmodeler:2.0.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.*")    , emit: db
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    BuildDatabase \\
        -name $prefix \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version | sed 's/RepeatModeler version //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nhr
    touch ${prefix}.nin
    touch ${prefix}.njs
    touch ${prefix}.nnd
    touch ${prefix}.nni
    touch ${prefix}.nog
    touch ${prefix}.nsq
    touch ${prefix}.translation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version | sed 's/RepeatModeler version //')
    END_VERSIONS
    """
}

process HELITRONSCANNER_DRAW {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/helitronscanner:1.0--hdfd78af_0':
        'biocontainers/helitronscanner:1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(head)
    tuple val(meta3), path(tail)

    output:
    tuple val(meta), path("*.draw")     , emit: draw
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: ''
    def args2       = task.ext.args2    ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"
    if ( !task.memory ) { error '[HELITRONSCANNER_DRAW] Available memory not known. Specify process memory requirements to fix this.' }
    def avail_mem   = (task.memory.giga*0.8).intValue()
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi

    HelitronScanner \\
        pairends \\
        -Xmx${avail_mem}g \\
        -head_score $head \\
        -tail_score $tail \\
        -output ${prefix}.pairends \\
        ${args2}

    HelitronScanner \\
        draw \\
        -Xmx${avail_mem}g \\
        -pscore ${prefix}.pairends \\
        -g $fasta \\
        -output ${prefix}.draw \\
        ${args}

    mv ${prefix}.draw.hel.fa \\
        ${prefix}.draw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        helitronscanner: \$(HelitronScanner |& sed -n 's/HelitronScanner V\\(.*\\)/V\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.draw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        helitronscanner: \$(HelitronScanner |& sed -n 's/HelitronScanner V\\(.*\\)/V\\1/p')
    END_VERSIONS
    """
}

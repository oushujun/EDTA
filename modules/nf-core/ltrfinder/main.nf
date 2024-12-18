process LTRFINDER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ltr_finder_parallel:1.1--hdfd78af_0':
        'biocontainers/ltr_finder_parallel:1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.scn")      , emit: scn
    tuple val(meta), path("*.gff3")     , emit: gff
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    LTR_FINDER_parallel \\
        -seq $fasta \\
        -threads $task.cpus \\
        $args

    mv "${fasta}.finder.combine.scn"    "${prefix}.scn"
    mv "${fasta}.finder.combine.gff3"   "${prefix}.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LTR_FINDER_parallel: \$(LTR_FINDER_parallel -h | grep 'Version:' | sed 's/Version: //')
        ltr_finder: \$(ltr_finder -h 2>&1 | grep 'ltr_finder' | sed 's/ltr_finder //')
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.scn"
    touch "${prefix}.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LTR_FINDER_parallel: \$(LTR_FINDER_parallel -h | grep 'Version:' | sed 's/Version: //')
        ltr_finder: \$(ltr_finder -h 2>&1 | grep 'ltr_finder' | sed 's/ltr_finder //')
    END_VERSIONS
    """
}

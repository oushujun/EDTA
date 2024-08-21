process ANNOSINE {
    tag "${meta.id}"
    label 'process_low' // maybe medium, previously had 12h and 16Gb ram
    
    input:
        tuple val(meta), path(assembly)

    output:
        tuple val(meta), path(assembly), path("${meta.id}.Seed_SINE.fa"), emit: seed_sine

        // Program does not output version, so must be hardcoded to what is specifed in conda/docker
        eval("echo '2.0.7'"), topic: versions
    
    conda 'bioconda::annosine2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annosine2%3A2.0.7--pyh7cba7a3_0':
        'quay.io/biocontainers/annosine2:2.0.7--pyh7cba7a3_0' }"
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch Seed_SINE.fa
    AnnoSINE_v2 -t ${task.cpus} \\
        -a 2 \\
        --num_alignments 50000 \\
        -rpm 0 \\
        --copy_number 3 \\
        --shift 100 \\
        $args \\
        -auto 1 3 ${assembly} ${prefix}

    mv Seed_SINE.fa ${prefix}.Seed_SINE.fa
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Seed_SINE.fa
    """
    
}
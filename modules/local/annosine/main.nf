process ANNOSINE {
    tag "${meta.id}"
    label 'process_low'
    
    input:
        tuple val(data), path(assembly)

    output:
        tuple val(data), path(assembly), path("${meta.id}.Seed_SINE.fa"), emit: seed_sine

        // todo test
        eval("AnnoSINE_v2 --version"), topic: versions
    
    conda 'bioconda::annosine2'
    container 'https://depot.galaxyproject.org/singularity/annosine2%3A2.0.7--pyh7cba7a3_0'

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
        -auto 1 3 ${data.assembly} .

    mv Seed_SINE.fa ${prefix}.Seed_SINE.fa
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Seed_SINE.fa
    """
    
}
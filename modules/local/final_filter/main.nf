process FINAL_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91c1946f5edb90aa9d2c26eb14e7348f8f4e94dda03bdb2185b72cb3c2d54932/data':
        'community.wave.seqera.io/library/blast_repeatmasker:816962ca420d1d16' }"

    input:
    tuple val(meta), path(genome, name: 'final/genome')
    path('final/genome.EDTA.raw.fa.cln')
    path('final/genome.EDTA.intact.fa.cln2')
    path('final/genome.EDTA.intact.raw.gff3')

    output:
    tuple val(meta), path('*.intact.fa')    , emit: intact_fa
    tuple val(meta), path('*.intact.gff3')  , emit: intact_gff
    tuple val(meta), path('*.TElib.fa')     , emit: telib_fa
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix          = task.ext.prefix   ?: "${meta.id}"
    """
    setup_final_filter.sh \\
        $task.cpus

    cd final
    
    perl final_filter.pl

    cd -

    mv \\
        final/genome.EDTA.intact.fa \\
        ${prefix}.intact.fa
    
    mv \\
        final/genome.EDTA.intact.gff3 \\
        ${prefix}.intact.gff3
    
    mv \\
        final/genome.EDTA.TElib.fa \\
        ${prefix}.TElib.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.intact.fa
    touch ${prefix}.intact.gff3
    touch ${prefix}.TElib.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """
}
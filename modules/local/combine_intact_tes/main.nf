process COMBINE_INTACT_TES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb105ba7d586ab31c1f39feb04c0255a39cc5a55ae7f6ea53f4bf76cdba8a3e5/data':
        'community.wave.seqera.io/library/mdust_tesorter_trf_perl:3424609103d3b065' }"

    input:
    tuple val(meta), path(genome, name: 'raw/genome')
    path "raw/genome.LTR.intact.raw.fa"
    path "raw/genome.LTR.intact.raw.gff3"
    path "raw/genome.TIR.intact.raw.fa"
    path "raw/genome.TIR.intact.raw.bed"
    path "raw/genome.Helitron.intact.raw.fa"
    path "raw/genome.Helitron.intact.raw.bed"

    output:
    tuple val(meta), path('*.intact.fa')    , emit: intact_raw_fa
    tuple val(meta), path('*.intact.gff3')  , emit: intact_raw_gff
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    setup_combine_intact_TEs.sh \\
        $task.cpus

    cd raw
    
    perl combine_intact_TEs.pl
    
    cd -

    mv \\
        raw/genome.EDTA.intact.raw.fa \\
        ${prefix}.intact.fa
    
    mv \\
        raw/genome.EDTA.intact.raw.gff3 \\
        ${prefix}.intact.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.intact.fa
    touch ${prefix}.intact.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
    END_VERSIONS
    """
}
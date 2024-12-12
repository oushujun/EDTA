process TIR_LEARNER_POSTPROCESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb105ba7d586ab31c1f39feb04c0255a39cc5a55ae7f6ea53f4bf76cdba8a3e5/data':
        'community.wave.seqera.io/library/mdust_tesorter_trf_perl:3424609103d3b065' }"

    input:
    tuple val(meta), path(genome, name: 'input/genome')
    path "input/TIR-Learner-Result/TIR-Learner_FinalAnn.fa"
    path "input/TIR-Learner-Result/TIR-Learner_FinalAnn.gff3"

    output:
    tuple val(meta), path('*.TIR.intact.raw.fa')    , emit: intact_raw_fa
    tuple val(meta), path('*.TIR.intact.raw.gff3')  , emit: intact_raw_gff3
    tuple val(meta), path('*.TIR.intact.raw.bed')   , emit: intact_raw_bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def MDUST_VERSION = '2006.10.17' // WARN: Manually update when changing Bioconda assets
    if ( "$prefix" == 'genome' ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    setup_TIR-Learner_postprocess.sh \\
        $task.cpus

    cd input
    
    perl TIR-Learner_postprocess.pl
    
    cd -

    mv \\
        genome.TIR.intact.raw.fa \\
        ${prefix}.TIR.intact.raw.fa
    
    mv \\
        genome.TIR.intact.raw.gff3 \\
        ${prefix}.TIR.intact.raw.gff3

    mv \\
        genome.TIR.intact.raw.bed \\
        ${prefix}.TIR.intact.raw.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        trf: \$(trf -v |& sed -n 's|.*Version \\(.*\\)|\\1|p')
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        mdust: $MDUST_VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def MDUST_VERSION = '2006.10.17' // WARN: Manually update when changing Bioconda assets
    if ( "$prefix" == 'genome' ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.TIR.intact.raw.fa
    touch ${prefix}.TIR.intact.raw.gff3
    touch ${prefix}.TIR.intact.raw.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        trf: \$(trf -v |& sed -n 's|.*Version \\(.*\\)|\\1|p')
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        mdust: $MDUST_VERSION
    END_VERSIONS
    """
}
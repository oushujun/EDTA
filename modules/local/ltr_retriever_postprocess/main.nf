process LTR_RETRIEVER_POSTPROCESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64d26063bedc2efcba20750a408a21f50907986d0d9aee685b03d7d05d3fbd8b/data':
        'community.wave.seqera.io/library/mdust_tesorter_trf_perl:3424609103d3b065' }"

    input:
    tuple val(meta), path(genome, name: 'input/genome')
    path "input/genome.pass.list"
    path "input/genome.pass.list.gff3"
    path "input/genome.defalse"
    path "input/genome.LTRlib.fa"

    output:
    tuple val(meta), path('*.LTR.raw.fa')                   , emit: raw_fa
    tuple val(meta), path('*.LTR.intact.raw.fa')            , emit: intact_raw_fa
    tuple val(meta), path('*.LTR.intact.raw.gff3')          , emit: intact_raw_gff3
    tuple val(meta), path('*.LTR.intact.raw.fa.anno.list')  , emit: anno_list
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def MDUST_VERSION = '2006.10.17' // WARN: Manually update when changing Bioconda assets
    if ( "$prefix" == 'genome' ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    setup_LTR_retriever_postprocess.sh \\
        $task.cpus

    cd input
    
    perl LTR_retriever_postprocess.pl
    
    cd -

    mv \\
        genome.LTR.raw.fa \\
        ${prefix}.LTR.raw.fa
    
    mv \\
        genome.LTR.intact.raw.fa \\
        ${prefix}.genome.LTR.intact.raw.fa

    mv \\
        genome.LTR.intact.raw.gff3 \\
        ${prefix}.LTR.intact.raw.gff3

    mv \\
        genome.LTR.intact.raw.fa.anno.list \\
        ${prefix}.LTR.intact.raw.fa.anno.list

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        trf: \$(trf -v |& sed -n 's|.*Version \\(.*\\)|\\1|p')
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        mdust: $MDUST_VERSION
    END_VERSIONS
    """
}
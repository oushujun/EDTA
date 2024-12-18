process HELITRONSCANNER_POSTPROCESS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb105ba7d586ab31c1f39feb04c0255a39cc5a55ae7f6ea53f4bf76cdba8a3e5/data':
        'community.wave.seqera.io/library/mdust_tesorter_trf_perl:3424609103d3b065' }"

    input:
    tuple val(meta), path(genome, name: 'input/genome')
    path "input/genome.HelitronScanner.filtered.fa"
    path "input/genome.HelitronScanner.filtered.ext.fa"

    output:
    tuple val(meta), path('*.Helitron.intact.anno.list')        , emit: anno_list
    tuple val(meta), path('*.Helitron.intact.bed')              , emit: raw_bed
    tuple val(meta), path('*.Helitron.intact.gff3')             , emit: raw_gff
    tuple val(meta), path('*.Helitron.intact.fa')               , emit: raw_fa
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def MDUST_VERSION = '2006.10.17' // WARN: Manually update when changing Bioconda assets
    """
    setup_HelitronScanner_postprocess.sh \\
        $task.cpus

    cd input
    
    perl HelitronScanner_postprocess.pl
    
    cd -

    mv \\
        genome.Helitron.intact.raw.fa.anno.list \\
        ${prefix}.Helitron.intact.anno.list

    mv \\
        genome.Helitron.intact.raw.bed \\
        ${prefix}.Helitron.intact.bed

    mv \\
        genome.Helitron.intact.raw.gff3 \\
        ${prefix}.Helitron.intact.gff3

    mv \\
        genome.Helitron.intact.raw.fa \\
        ${prefix}.Helitron.intact.fa

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
    """
    touch ${prefix}.Helitron.intact.anno.list
    touch ${prefix}.Helitron.intact.bed
    touch ${prefix}.Helitron.intact.gff3
    touch ${prefix}.Helitron.intact.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        trf: \$(trf -v |& sed -n 's|.*Version \\(.*\\)|\\1|p')
        TEsorter: \$(TEsorter -v | tr -d 'TEsorter ')
        mdust: $MDUST_VERSION
    END_VERSIONS
    """
}
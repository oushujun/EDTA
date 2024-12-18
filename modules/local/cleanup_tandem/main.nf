process CLEANUP_TANDEM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64d26063bedc2efcba20750a408a21f50907986d0d9aee685b03d7d05d3fbd8b/data':
        'community.wave.seqera.io/library/trf_perl:44f9844b1bf765c3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta")    , emit: cleaned_fasta
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args ?: ''
    if ( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    cleanup_tandem.pl \\
        -f $fasta \\
        $args \\
        > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        trf: \$(trf -v |& sed -n 's|.*Version \\(.*\\)|\\1|p')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    if ( "$fasta" == "${prefix}.fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        trf: \$(trf -v |& sed -n 's|.*Version \\(.*\\)|\\1|p')
    END_VERSIONS
    """
}
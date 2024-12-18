process FORMAT_HELITRONSCANNER_OUT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c0905ae7aa2a1c4135b448b49dd205b4eb370d4733b832594563ee06832ebd09/data':
        'community.wave.seqera.io/library/perl:67db1a8d53d64b14' }"

    input:
    tuple val(meta), path(genome)
    path hel_fa
    path rc_hel_fa

    output:
    tuple val(meta), path("*.HelitronScanner.filtered.tabout")  , emit: filtered_tabout
    tuple val(meta), path("*.HelitronScanner.filtered.fa")      , emit: filtered_fa     , optional: true
    tuple val(meta), path("*.HelitronScanner.filtered.ext.fa")  , emit: filtered_ext_fa , optional: true
    tuple val(meta), path("*.format_helitronscanner_out.log")   , emit: log
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "$genome"
    def args    = task.ext.args ?: ''
    """
    # Create symbolic links to match expected filenames
    ln -s $hel_fa ${genome}.HelitronScanner.draw.hel.fa
    ln -s $rc_hel_fa ${genome}.HelitronScanner.draw.rc.hel.fa

    format_helitronscanner_out.pl \\
        -genome $genome \\
        $args \\
        &> >(tee "${prefix}.format_helitronscanner_out.log" 2>&1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "$genome"
    """
    touch ${prefix}.HelitronScanner.filtered.tabout
    touch ${prefix}.format_helitronscanner_out.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
    END_VERSIONS
    """
}

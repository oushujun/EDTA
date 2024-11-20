process FORMAT_HELITRONSCANNER_OUT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-bioperl:1.7.8--hdfd78af_1':
        'biocontainers/perl-bioperl:1.7.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(genome)
    path hel_fa       // HelitronScanner.draw.hel.fa file
    path rc_hel_fa    // HelitronScanner.draw.rc.hel.fa file

    output:
    tuple val(meta), path("${meta.id}.HelitronScanner.filtered.tabout"), emit: filtered_tabout
    tuple val(meta), path("${meta.id}.HelitronScanner.filtered.fa"),     emit: filtered_fa, optional: true
    tuple val(meta), path("${meta.id}.HelitronScanner.filtered.ext.fa"), emit: filtered_ext_fa, optional: true
    tuple val(meta), path("${meta.id}.format_helitronscanner_out.log"),  emit: log
    path "versions.yml",                                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Set default values for the optional arguments
    def sitefilter   = task.ext.sitefilter   != null ? task.ext.sitefilter   : 1
    def minscore     = task.ext.minscore     != null ? task.ext.minscore     : 12
    def keepshorter  = task.ext.keepshorter  != null ? task.ext.keepshorter  : 1
    def extlen       = task.ext.extlen       != null ? task.ext.extlen       : 30
    def extout       = task.ext.extout       != null ? task.ext.extout       : 1
    def prefix       = meta.id

    """
    # Create symbolic links to match expected filenames
    ln -s $hel_fa ${genome}.HelitronScanner.draw.hel.fa
    ln -s $rc_hel_fa ${genome}.HelitronScanner.draw.rc.hel.fa

    # Run the Perl script with the provided arguments
    perl format_helitronscanner_out.pl \\
        -genome $genome \\
        -sitefilter $sitefilter \\
        -minscore $minscore \\
        -keepshorter $keepshorter \\
        -extlen $extlen \\
        -extout $extout \\
        &> >(tee "${prefix}.format_helitronscanner_out.log" 2>&1)

    # Capture version information
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl: \$(perl -v | grep 'This is perl' | awk '{print \$4}')
    END_VERSIONS
    """
}

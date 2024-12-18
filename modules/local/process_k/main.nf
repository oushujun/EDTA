process PROCESS_K {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/91c1946f5edb90aa9d2c26eb14e7348f8f4e94dda03bdb2185b72cb3c2d54932/data':
        'community.wave.seqera.io/library/blast_repeatmasker:816962ca420d1d16' }"

    input:
    tuple val(meta)     , path(genome, name: 'genome')
    path(ltr            , name: "genome.EDTA.raw/genome.LTR.raw.fa")
    path(ltrint         , name: "genome.EDTA.raw/genome.LTR.intact.raw.fa")
    path(sine           , name: "genome.EDTA.raw/genome.SINE.raw.fa")
    path(line           , name: "genome.EDTA.raw/genome.LINE.raw.fa")
    path(tir            , name: "genome.EDTA.raw/genome.TIR.intact.raw.fa")
    path(helitron       , name: "genome.EDTA.raw/genome.Helitron.intact.raw.fa")

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: stg1_fa
    tuple val(meta), path('*.intact.fasta') , emit: intact_fa
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix              = task.ext.prefix   ?: "${meta.id}"
    def args            = task.ext.args     ?: ''

    def touch_ltr       = ltr               ? ''    : "touch genome.EDTA.raw/genome.LTR.raw.fa"
    def touch_ltrint    = ltrint            ? ''    : "touch genome.EDTA.raw/genome.LTR.intact.raw.fa"
    def touch_sine      = sine              ? ''    : "touch genome.EDTA.raw/genome.SINE.raw.fa"
    def touch_line      = line              ? ''    : "touch genome.EDTA.raw/genome.LINE.raw.fa"
    def touch_tir       = tir               ? ''    : "touch genome.EDTA.raw/genome.TIR.intact.raw.fa"
    def touch_helitron  = helitron          ? ''    : "touch genome.EDTA.raw/genome.Helitron.intact.raw.fa"
    """
    mkdir -p genome.EDTA.raw
    $touch_ltr
    $touch_ltrint
    $touch_sine
    $touch_line
    $touch_tir
    $touch_helitron

    PROGRAM_PATH=\$(dirname \$(which EDTA_processK.pl))

    sed "s|\\\$script_path/bin|\$PROGRAM_PATH|g" \\
        \$PROGRAM_PATH/EDTA_processK.pl \\
        > EDTA_processK.pl
    
    perl \\
        EDTA_processK.pl \\
        -genome genome \\
        -ltr genome.EDTA.raw/genome.LTR.raw.fa \\
        -ltrint genome.EDTA.raw/genome.LTR.intact.raw.fa \\
        -sine genome.EDTA.raw/genome.SINE.raw.fa \\
        -line genome.EDTA.raw/genome.LINE.raw.fa \\
        -tir genome.EDTA.raw/genome.TIR.intact.raw.fa \\
        -helitron genome.EDTA.raw/genome.Helitron.intact.raw.fa \\
        $args \\
        ${task.cpus}
    
    mv \\
        genome.EDTA.combine/genome.EDTA.fa.stg1 \\
        ${prefix}.fasta
    
    mv \\
        genome.EDTA.combine/genome.EDTA.intact.fa.cln \\
        ${prefix}.intact.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """

    stub:
    prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    touch ${prefix}.intact.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -v | sed -n 's|This is perl.*(\\(.*\\)).*|\\1|p')
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """
}
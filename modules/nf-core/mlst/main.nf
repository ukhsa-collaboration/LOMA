process MLST {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/mlst:2.23.0--hdfd78af_0' :
        'https://depot.galaxyproject.org/singularity/mlst:2.23.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    scheme = meta.mlst_scheme != "" ? "--scheme ${meta.mlst_scheme}" : ""
    additional_params = (meta.mlst_scheme.matches("yersinia_mcnally") && params.MLST.yersinia_blastdb && params.MLST.yersinia_datadir) ? "--blastdb ${params.MLST.yersinia_blastdb} --datadir ${params.MLST.yersinia_datadir}" : ""

    """
    mlst \\
        $args \\
        $scheme \\
        $additional_params \\
        --threads $task.cpus \\
        $fasta \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
    """

}

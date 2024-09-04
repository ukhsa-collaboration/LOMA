process METADECODER_SEED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/metadecoder:latest' :
        'quay.io/djberger/metadecoder:latest' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.seed")   , emit: seed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metadecoder \\
       seed \\
       --threads $task.cpus \\
       -f $fasta \\
       -o ${prefix}.metadecoder.seed 2>/dev/null || echo "Metadecoder failed, error is ignored"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        METADECODER: \$( metadecoder --version | sed 's/metadecoder //g' )
    END_VERSIONS
"""
}

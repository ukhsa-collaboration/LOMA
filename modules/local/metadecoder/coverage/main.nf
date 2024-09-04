process METADECODER_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/metadecoder:latest' :
        'quay.io/djberger/metadecoder:latest' }"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.coverage")   , emit: cov
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metadecoder \\
       coverage \\
       -s $sam \\
       -o ${prefix}.metadecoder.coverage 2>/dev/null || echo "Metadecoder failed, error is ignored"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaDecoder: \$( metadecoder --version | sed 's/metadecoder //g' )
    END_VERSIONS
"""
}

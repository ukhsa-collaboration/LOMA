process METADECODER_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/metadecoder:latest' :
        'quay.io/djberger/metadecoder:latest' }"

    input:
    tuple val(meta), path(fasta), path(cov), path(seed)

    output:
    tuple val(meta), path("*_output_bins/*.fasta")   , emit: binned_fastas
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_output_bins
    touch ${prefix}_output_bins/empty.fasta

    metadecoder \\
       cluster \\
       --min_sequence_length 2000 \\
       $args \\
       -f $fasta \\
       -c ${prefix}.metadecoder.coverage \\
       -s ${prefix}.metadecoder.seed \\
       -o ${prefix}_output_bins/${prefix} 2>/dev/null || echo "Metadecoder failed, error is ignored"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaDecoder: \$( metadecoder --version | sed 's/metadecoder //g' )
    END_VERSIONS
"""
}

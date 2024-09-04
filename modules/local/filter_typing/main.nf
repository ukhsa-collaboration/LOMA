process FILTER_TYPING {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(seqs), path(results1), path(results2)
    val(mode)

    output:
    tuple val(meta), path(seqs), path("*.ftype.tsv"), emit: typing
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = mode.matches("ecoli") ? "--mlst $results1 --krocus $results2" : "--seqsero2 $results1 --sistr $results2"

    """
    filter_typing.py \\
       $input \\
       --output $prefix \\
       --bin ${meta.bin}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_typing.py: \$(filter_typing.py --version 2>&1 | tail -1 | cut -f2 -d " ")
    END_VERSIONS

    """
}

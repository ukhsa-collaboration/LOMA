process SEQSERO2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/seqsero2:1.3.1--pyhdfd78af_1' :
        'https://depot.galaxyproject.org/singularity/seqsero2:1.3.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(seqs)

    output:
    tuple val(meta), path("results/*_log.txt")   , emit: log
    tuple val(meta), path("*.result.tsv"), emit: tsv
    tuple val(meta), path("*.result.txt"), emit: txt
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    SeqSero2_package.py \\
        $args \\
        -d results/ \\
        -n $prefix \\
        -p $task.cpus \\
        -i $seqs

    mv results/*_result.tsv ${prefix}.result.tsv
    mv results/*_result.txt ${prefix}.result.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqsero2: \$( echo \$( SeqSero2_package.py --version 2>&1) | sed 's/^.*SeqSero2_package.py //' )
    END_VERSIONS
    """
}

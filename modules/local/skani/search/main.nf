process SKANI_SEARCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/skani:0.2.1--h4ac6f70_0' :
        'quay.io/biocontainers/skani:0.2.1--h4ac6f70_0' }"

    input:
    tuple val(meta), path(fasta)
    path("database")

    output:
    tuple val(meta), path("${prefix}.skani_results.txt")         , emit: summary
    path("versions.yml")                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    skani \\
      search \\
      $fasta \\
      --qi \\
      -d ${database} \\
      $args \\
      -o ${prefix}.skani_results.txt \\
      -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        skani: \$(skani -V | sed "s/skani //")
    END_VERSIONS
    """
}

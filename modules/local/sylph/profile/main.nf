process SYLPH_PROFILE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/sylph:0.6.1--h4ac6f70_0':
        'biocontainers/sylph:0.6.1--h4ac6f70_0' }"

    input:
    tuple val(meta), path(sketch)
    path(db)

    output:
    tuple val(meta), path('*.tsv'), optional: true, emit: results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    sylph profile \\
       $args \\
       -u \\
       -t $task.cpus \\
       -o ${prefix}.sylph_results.tsv \\
       ${sketch} \\
       $db

    find . -type f -empty -print -delete

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(echo \$( sylph -V 2>&1 | sed 's/sylph //' ))
    END_VERSIONS
    """
}

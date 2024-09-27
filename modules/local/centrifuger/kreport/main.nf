process CENTRIFUGER_KREPORT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/centrifuger:1.0.3--hdcf5f25_0' :
        'https://depot.galaxyproject.org/singularity/centrifuger:1.0.3--hdcf5f25_0' }"

    input:
    tuple val(meta), path(results)
    path(db)

    output:
    tuple val(meta), path('*.kreport.txt')                 , emit: kreport
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    db_name=`find -L ${db} -name "*.1.cfr" -not -name "._*"  | sed 's/\\.1.cfr\$//'`

    centrifuger-kreport \\
        -x \$db_name \\
        $results \\
        $args > ${prefix}.kreport.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuger: \$(  centrifuger -v  | sed -n 1p | sed 's/Centrifuger v//')
    END_VERSIONS
    """
}

process SHIGATYPER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/shigatyper:2.0.5--pyhdfd78af_0' :
        'https://depot.galaxyproject.org/singularity/shigatyper:2.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.summary.tsv")     , optional: true, emit: tsv
    tuple val(meta), path("*-hits.tsv"), optional: true, emit: hits
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.is_ont) {
        """
        shigatyper \\
            $args \\
            --SE $reads \\
            --ont \\
            --name $prefix 2>/dev/null || echo ""

        mv *-hits.tsv ${prefix}-hits.txt
        mv *.tsv ${prefix}.summary.tsv 
        mv ${prefix}-hits.txt ${prefix}-hits.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    } else if (meta.single_end) {
        """
        shigatyper \\
            $args \\
            --SE $reads \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    } else {
        """
        shigatyper \\
            $args \\
            --R1 ${reads[0]} \\
            --R2 ${reads[1]} \\
            --name $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shigatyper: \$(echo \$(shigatyper --version 2>&1) | sed 's/^.*ShigaTyper //' )
        END_VERSIONS
        """
    }
}

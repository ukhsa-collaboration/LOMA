process FILTER_BLAST {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(workflow), path(blast_results)

    output:
    tuple val(meta), path("*_filtered_hits.tsv"), emit: filtered_hits
    tuple val(meta), path("*_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_blasthits.py \\
       --workflow workflow.txt \\
       --blast_results ${blast_results} \\
       --output $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_blasthits.py: \$(filter_blasthits.py --version 2>&1 | tail -1 | cut -f2 -d " ")
    END_VERSIONS

    """
}

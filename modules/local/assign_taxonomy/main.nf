process ASSIGN_TAXONOMY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path("input_bins/*"), path(gtdb_results)
    tuple val(gtdb_ani_cutoff), val(gtdb_aln_frac)
    path(gtdb_definition_table)

    output:
    tuple val(meta), path("${meta.id}.gtdb_filtered_report.txt"), optional: true, emit: report
    tuple val(meta), path("*.fasta"), emit: all_fasta, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_gtdbtk.py \\
       --input *.gtdbtk.bac120.* \\
       --output ${prefix} \\
       --ANI_cutoff $gtdb_ani_cutoff \\
       --align_fraction_cutoff $gtdb_aln_frac \\
       --mode Nanopore \\
       --gtdb_definition_table $gtdb_definition_table 2>/dev/null || echo ""

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_gtdbtk.py.py: \$(filter_gtdbtk.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS
    """
}

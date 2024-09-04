process SUMMARIZE_AMR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(input_files)
    path(template)

    output:
    tuple val(meta), path("*.html"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_1 = input_files.toString().contains(".summarize.hamronization_combined_report.tsv") ? "--hamronization_summary *.summarize.hamronization_combined_report.tsv" : ""
    def in_2 = input_files.toString().contains(".bin_summary.tsv") ? "--bin_summary *.bin_summary.tsv" : ""

    """
    write_amr_report.py \\
       $in_1 \\
       $in_2 \\
       --report_template $template \\
       --logo $params.logo \\
       --sample_id $meta.id \\
       --run_id $meta.run_id \\
       --barcode $meta.barcode \\
       --sample_type $meta.target \\
       --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_amr_report.py: \$(write_amr_report.py --version 2>&1 | tail -1 | cut -f2 -d " ")
    END_VERSIONS
    """
}

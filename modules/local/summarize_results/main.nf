process SUMMARIZE_RESULTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(input_files, stageAs: "RESULTS/*")
    path(template)

    output:
    path("*.summary_report.html"), emit: summary_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def in_1 = input_files.toString().contains(".host_reads.txt") ? "--hostreads RESULTS/*.host_reads.txt" : ""
    def in_2 = input_files.toString().contains(".nanostats.tsv") ? "--readqc RESULTS/*.nanostats.tsv" : ""
    def in_3 = input_files.toString().contains(".bracken_kraken2.target_species.tsv") ? "--kraken2_bracken_targets RESULTS/*.bracken_kraken2.target_species.tsv" : ""
    def in_4 = input_files.toString().contains(".bracken_kraken2.tsv") ? "--kraken2_bracken_all RESULTS/*.bracken_kraken2.tsv" : ""
    def in_5 = input_files.toString().contains(".contig_summary.tsv") ? "--contig_summary RESULTS/*.contig_summary.tsv" : ""
    def in_6 = input_files.toString().contains(".bin_summary.tsv") ? "--bin_summary RESULTS/*.bin_summary.tsv" : ""
    def in_7 = input_files.toString().contains(".krocus_cc.tsv") ? "--krocus RESULTS/*.krocus_cc.tsv" : ""
    def in_8 = input_files.toString().contains(".mlst_cc.tsv") ? "--mlst RESULTS/*.mlst_cc.tsv" : ""
    def in_9 = input_files.toString().contains(".gene_hits.summary.tsv") ? "--genefinding RESULTS/*.gene_hits.summary.tsv" : ""
    def in_10 = input_files.toString().contains(".plasmidfinder.") ? "--plasmidfinder RESULTS/*.plasmidfinder.*" : ""
    def in_11 = input_files.toString().contains(".summarize.hamronization_combined_report.tsv") ? "--hamronization_summary RESULTS/*.summarize.hamronization_combined_report.tsv" : ""
    def in_12 = input_files.toString().contains(".salmonella_typing_short.tsv") ? "--salmonella_typing RESULTS/*.salmonella_typing_short.tsv" : ""
    def in_13 = input_files.toString().contains("lissero") ? "--lmonocytogenes_typing RESULTS/*.lissero.tsv" : ""
    def in_14 = input_files.toString().contains(".ecoli_typing_short.tsv") ? "--ecoli_typing RESULTS/*.ecoli_typing_short.tsv" : ""
    def in_15 = input_files.toString().contains(".bracken_centrifuger.target_species.tsv") ? "--centrifuger_bracken_targets RESULTS/*.bracken_centrifuger.target_species.tsv" : ""
    def in_16 = input_files.toString().contains(".bracken_centrifuger.tsv") ? "--centrifuger_bracken_all RESULTS/*.bracken_centrifuger.tsv" : ""


    """

    write_summary_report.py \\
       --output $prefix \\
       --sample_id $meta.id \\
       --run_id $meta.run_id \\
       --barcode $meta.barcode \\
       --sample_type $meta.target \\
       --logo $params.logo \\
       --pipeline_version $params.pipeline_version \\
       $in_1 \\
       $in_2 \\
       $in_3 \\
       $in_4 \\
       $in_5 \\
       $in_6 \\
       --min_read_count $params.PARSE_KRAKEN2HITS.min_target_reads \\
       --min_read_prop $params.PARSE_KRAKEN2HITS.min_target_fraction \\
       $in_7 \\
       $in_8 \\
       $in_9 \\
       $in_10 \\
       $in_11 \\
       --amr_tool $params.CREATE_REPORT.amr_tool \\
       $in_12 \\
       $in_13 \\
       $in_14 \\
       $in_15 \\
       $in_16 \\
       --tax_mode $params.SUMMARIZE_RESULTS.tax_mode \\
       --report_template $template

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_summary_report.py: \$(write_summary_report.py --version 2>&1 | tail -1 | cut -f2 -d " ")
    END_VERSIONS
    """
}

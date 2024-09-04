process PLOT_TAXHITS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(taxpasta_kraken), path(taxpasta_centrifuger), path(sylph), path(kraken2_bracken), path(centrifuger_bracken)
    path(syl_fn)
    path(template)

    output:
    tuple val(meta), path("*.html"), emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_taxhits.py \\
       --taxpasta_kraken2 $taxpasta_kraken \\
       --taxpasta_centrifuger $taxpasta_centrifuger \\
       --sylph $sylph \\
       --bracken_kraken2 $kraken2_bracken \\
       --bracken_centrifuger $centrifuger_bracken \\
       --syl_fn $syl_fn \\
       --logo $params.logo \\
       --sample_id ${meta.id} \\
       --run_id ${meta.run_id} \\
       --barcode ${meta.barcode} \\
       --sample_type ${meta.target} \\
       --report_template $template \\
       --output $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plot_taxhits.py: \$(plot_taxhits.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS

    """
}

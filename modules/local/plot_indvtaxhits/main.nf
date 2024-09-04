process PLOT_INDVTAXHITS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(taxhits)
    val(mode)
    val(db)
    path(template)

    output:
    tuple val(meta), path("*.html"), emit: report
    path "versions.yml", emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxp_kraken2 = mode.matches("Kraken2, Taxpasta") ? "--taxpasta_kraken2 $taxhits" : ""
    def taxp_centrifuger = mode.matches("Centrifuger, Taxpasta") ? "--taxpasta_centrifuger $taxhits" : ""
    def bracken_kraken2 = mode.matches("Kraken2, Bracken") ? "--bracken_kraken2 $taxhits" : ""
    def bracken_centrifuger = mode.matches("Centrifuger, Bracken") ? "--bracken_centrifuger $taxhits" : ""
    def sylph = mode.matches("Sylph") ? "--sylph $taxhits --syl_fn $db" : ""

    """
    plot_taxhits.py \\
       $args \\
       $taxp_kraken2 \\
       $taxp_centrifuger \\
       $bracken_kraken2 \\
       $bracken_centrifuger \\
       $sylph \\
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

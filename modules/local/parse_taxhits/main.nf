process PARSE_TAXHITS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(bracken_taxhits)
    path(targets)
    tuple val(min_reads), val(min_frac)
    val(mode)

    output:
    tuple val(meta), path("*.target_species.tsv"), emit: targets_filtered
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def opt = mode.matches("Sylph") ? "--gtdb_fn $params.TAXONOMIC_PROFILING.gtdb_metadata" : ""

    """
    parse_taxonomic_hits.py \\
       $opt \\
       --taxhits $bracken_taxhits \\
       --output ${prefix} \\
       --targets $targets \\
       --min_reads $min_reads \\
       --min_frac $min_frac \\
       --mode $mode

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_taxonomic_hits.py: \$( parse_taxonomic_hits.py --version | cut -f2 -d " " )
    END_VERSIONS

    """
}

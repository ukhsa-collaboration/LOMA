process FILTER_ECOLI {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(shigeifinder), path(shigatyper), path(ectyper), path(stecfinder), path(mlst), path(krocus), path(mykrobe)

    output:
    tuple val(meta), path("*_typing_*.tsv"), emit: typing
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def opt = mykrobe.toString().matches("") ? "" : "--mykrobe $mykrobe"

    """
    filter_ecoli.py \\
       --bin ${meta.bin} \\
       --shigeifinder $shigeifinder \\
       --shigatyper $shigatyper \\
       --ectyper $ectyper \\
       --stecfinder $stecfinder \\
       $opt \\
       --mlst $mlst \\
       --krocus $krocus \\
       --output ${prefix} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_ecoli.py: \$(filter_ecoli.py --version 2>&1 | tail -1 | cut -f2 -d " ")
    END_VERSIONS

    """
}

process CONVERT_MLSTTOCC {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(mlst)
    val(profiler)
    path(gastro_cc_db)

    output:
    tuple val(meta), path("*_cc.tsv"), emit: mlstwithcc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def krocus = profiler.matches("mlst") ? "" : "--krocus_scheme '${meta.krocus_scheme}'"

    """
    convert_mlsttocc.py \\
       --profile $mlst \\
       --clonal_complexes $gastro_cc_db \\
       --output ${meta.id} \\
       --sample ${meta.id} \\
       --bin ${meta.bin} \\
       --mode $profiler $krocus


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert_mlsttocc.py: \$(rename_bins.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS
    """
}

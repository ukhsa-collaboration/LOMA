process FETCH_FASTQ {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("*.fastq.gz"), optional: true, emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    barcode_in = (meta.barcode =~ /^RB/ || "$meta.barcode" =~ /^NB/) ? "--barcode ${meta.barcode}" : ""

    """
    process_basecalled.py \\
       --input-dir $params.run_dir \\
       --run_id ${meta.run_id} \\
       --output_prefix ${prefix} \\
       --folder_name guppy ${barcode_in}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        process_basecalled.py: \$(process_basecalled.py --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS
    """
}

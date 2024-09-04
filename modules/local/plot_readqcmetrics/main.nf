process PLOT_READQCMETRICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(preqc_fw), path(preqc_rv), path(preqc_nanostats), path(preqc_nanoplot_raw), path(postqc_fw), path(postqc_rv), path(postqc_nanostats), path(postqc_nanoplot_raw)
    path(template)

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*.nanostats.tsv"), emit: nanostats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_readqc.py \\
       --nanoplot_raw_pre $preqc_nanoplot_raw \\
       --nanoplot_raw_post $postqc_nanoplot_raw \\
       --nanostats_pre $preqc_nanostats \\
       --nanostats_post $postqc_nanostats \\
       --nucl_comp_pre_fw $preqc_fw \\
       --nucl_comp_pre_rv $preqc_rv \\
       --nucl_comp_post_fw $postqc_fw \\
       --nucl_comp_post_rv $postqc_rv \\
       --report_template $template \\
       --logo $params.logo \\
       --sample_id $meta.id \\
       --run_id $meta.run_id \\
       --barcode $meta.barcode \\
       --sample_type $meta.target \\
       --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plot_readqc.py: \$(plot_readqc.py --version 2>&1 | tail -1 | cut -f2 -d " ")
    END_VERSIONS
    """
}

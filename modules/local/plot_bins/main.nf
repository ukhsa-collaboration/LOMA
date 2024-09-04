process PLOT_BINS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/djberger/lma_img:latest' :
        'quay.io/djberger/lma_img:latest' }"

    input:
    tuple val(meta), path(cov), path(fstats), path(skani), path(plasmid_score), path(checkm), path(bin_tax), path(asm_stats)
    path(syl_fn)
    path(template)

    output:
    tuple val(meta), path("*.html"), path("*.html"), emit: html 
    tuple val(meta), path("*.contig_summary.tsv"), emit: contig_summary
    tuple val(meta), path("*.bin_summary.tsv"), emit: bin_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_bins.py \\
       $args \\
       --asm_stats $asm_stats \\
       --cov $cov \\
       --fstats $fstats \\
       --bin_tax *.bac*.summary.tsv \\
       --checkm $checkm \\
       --skani $skani \\
       --genomad_plasmid $plasmid_score \\
       --output $prefix \\
       --report_template $template \\
       --logo $params.logo \\
       --sample_id $meta.id \\
       --run_id $meta.run_id \\
       --barcode $meta.barcode \\
       --sample_type $meta.target \\
       --gtdb_fn $syl_fn 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plot_bins.py: \$(plot_bins.py --version 2>&1 | tail -1 | cut -f2 -d " ")
    END_VERSIONS

    """
}

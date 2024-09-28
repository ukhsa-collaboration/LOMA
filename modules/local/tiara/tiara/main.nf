process TIARA_TIARA {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/tiara:1.0.3' :
        'https://depot.galaxyproject.org/singularity/tiara:1.0.3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.{txt,txt.gz}")  , emit: classifications
    tuple val(meta), path("log_*.{txt,txt.gz}")      , emit: log
    tuple val(meta), path("*.{fasta,fasta.gz}")          , emit: fasta, optional: true
    tuple val(meta), path("*.tiara_headers.txt"), emit: prokaryote_tiara, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    tiara -i ${fasta} \
        -o ${prefix}.txt \
        --threads ${task.cpus} \
        ${args}

    awk -F'\t' '\$2 == "bacteria" || \$2 == "archaea" || \$2 == "unknown" || \$2=="prokarya"' OFS='\t' ${prefix}.txt | cut -f1 > ${prefix}.tiara_headers.txt

    ## fix gzip flag weirdness and ensure consistent .fasta filename output
    ## check if fasta files are being output
    if echo "${args}" | grep -qE "tf|to-fasta"; then
        ## check if we've asked for gzip output, then rename files consistently
        if echo "${args}" | grep -q "gz"; then
            find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}.fasta.gz"' \\;
        else
            find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}.fasta"' \\;
        fi
    fi


    find . -type f -empty -print -delete

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.out.txt
    touch log_${prefix}.out.txt
    touch bacteria_${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """
}

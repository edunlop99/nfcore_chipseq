process HOMER_FINDMOTIFS {
    tag "${meta.id}"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(peaks)
    path fasta

    output:
    tuple val(meta), path("${prefix}"),                   emit: motifs
    tuple val(meta), path("${prefix}/knownResults.txt"),  emit: known
    tuple val("${task.process}"), val('homer'), val("4.11"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_motifs"
    """
    mkdir -p preparsed

    findMotifsGenome.pl \\
        ${peaks} \\
        ${fasta} \\
        ${prefix} \\
        -p ${task.cpus} \\
        -preparsedDir preparsed \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_motifs"
    """
    mkdir -p ${prefix}
    touch ${prefix}/knownResults.txt
    """
}
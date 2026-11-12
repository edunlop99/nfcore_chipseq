// DiffBind for a single contrast: two groups of one antibody.
process DIFFBIND {
    tag "${meta.id}"
    label 'process_high'

    // TODO: pin a DiffBind image. Tags: https://quay.io/repository/biocontainers/bioconductor-diffbind?tab=tags
    conda "bioconda::bioconductor-diffbind"
    container "quay.io/biocontainers/bioconductor-diffbind:<version--build>"

    input:
    tuple val(meta), val(metas), path(peaks), path(ip_bams), path(ip_bais),
          val(ctrl_names), path(ctrl_bams), path(ctrl_bais)                 // control lists are [] without controls
    val is_narrow_peak

    output:
    tuple val(meta), path("*.samplesheet.csv")      , emit: samplesheet
    tuple val(meta), path("*.results.csv")          , emit: results , optional: true
    tuple val(meta), path("*.consensus_scores.csv") , emit: scores  , optional: true
    tuple val(meta), path("*.plots.pdf")            , emit: plots   , optional: true
    tuple val(meta), path("*.dba.rds")              , emit: rds     , optional: true
    tuple val(meta), path("*.log.txt")              , emit: log
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def has_ctrl = ctrl_names.any()

    // DiffBind sample sheet; control columns only when there are controls
    def factor_cols = [ Tissue: 'tissue', Condition: 'condition', Treatment: 'treatment' ].findAll { col, key -> metas.any { it[key] } } // factor columns that are blank for every sample are left out of the sheet
    def cols = [ 'SampleID' ] + factor_cols.keySet() + [ 'Replicate' ] +
            ( has_ctrl ? [ 'ControlID', 'bamControl' ] : [] ) +
            [ 'bamReads', 'Peaks' ]
    def rows = metas.withIndex().collect { m, i ->
        ( [ m.id ] + factor_cols.values().collect { m[it] } + [ m.replicate ] +
        ( has_ctrl ? [ m.control, ctrl_names[i] ] : [] ) +
        [ ip_bams[i].name, peaks[i].name ] ).join(',')
    }
    // One quoted argument per line; printf repeats '%s\n' for each
    def quoted = ([ cols.join(',') ] + rows).collect { "'${it}'" }.join(' ')
    """
    set -o pipefail

    printf '%s\\n' ${quoted} > ${prefix}.samplesheet.csv

    diffbind.R \\
        --samplesheet=${prefix}.samplesheet.csv \\
        --prefix=${prefix} \\
        --numerator=${meta.numerator} \\
        --denominator=${meta.denominator} \\
        --peak_type=${is_narrow_peak ? 'narrow' : 'broad'} \\
        --cores=${task.cpus} \\
        ${args} \\
        2>&1 | tee ${prefix}.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
        bioconductor-diffbind: \$(Rscript -e 'suppressMessages(library(DiffBind)); cat(as.character(packageVersion("DiffBind")))')
    END_VERSIONS
    """

    stub:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def has_ctrl = ctrl_names.any()
    def cols = [ 'SampleID', 'Tissue', 'Condition', 'Treatment', 'Replicate' ] +
               ( has_ctrl ? [ 'ControlID', 'bamControl' ] : [] ) + [ 'bamReads', 'Peaks' ]
    def rows = metas.withIndex().collect { m, i ->
        ( [ m.id, m.tissue, m.condition, m.treatment, m.replicate ] +
          ( has_ctrl ? [ m.control, ctrl_names[i] ] : [] ) + [ ip_bams[i].name, peaks[i].name ] ).join(',')
    }
    def quoted = ([ cols.join(',') ] + rows).collect { "'${it}'" }.join(' ')
    """
    printf '%s\\n' ${quoted} > ${prefix}.samplesheet.csv
    touch ${prefix}.results.csv ${prefix}.consensus_scores.csv ${prefix}.plots.pdf ${prefix}.dba.rds ${prefix}.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-diffbind: stub
    END_VERSIONS
    """
}

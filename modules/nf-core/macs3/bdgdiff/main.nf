// MACS3 bdgdiff between two groups' pooled bedGraphs.
// Outputs: <prefix>_cond1.bed (higher in cond1), <prefix>_cond2.bed (higher in cond2),
//          <prefix>_common.bed, and <prefix>.contrast.tsv recording which group is cond1/cond2.
process MACS3_BDGDIFF {
    tag "${meta.id}"
    label 'process_medium'

    // TODO: pin the same MACS3 version the main workflow uses
    conda "bioconda::macs3"
    container "quay.io/biocontainers/macs3:<version--build>"

    input:
    tuple val(meta), path(t1), path(c1), path(xls1), path(t2), path(c2), path(xls2)

    output:
    tuple val(meta), path("*_cond1.bed")    , emit: cond1
    tuple val(meta), path("*_cond2.bed")    , emit: cond2
    tuple val(meta), path("*_common.bed")   , emit: common
    tuple val(meta), path("*.contrast.tsv") , emit: contrast
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Depth (millions of reads) = smaller of treatment/control after filtering,
    # because callpeak scales the larger library down to the smaller one by default.
    # Paired-end (BAMPE) runs report "fragments" instead of "tags".
    depth() {
        local t c
        t=\$(grep -E '^# (tags|fragments) after filtering in treatment' "\$1" | awk -F': ' '{print \$2}')
        c=\$(grep -E '^# (tags|fragments) after filtering in control'   "\$1" | awk -F': ' '{print \$2}')
        if [ -z "\$t" ]; then echo "Could not read depth from \$1" >&2; exit 1; fi
        [ -z "\$c" ] && c=\$t
        awk -v t="\$t" -v c="\$c" 'BEGIN { m = (t < c ? t : c); printf "%.6f", m / 1000000 }'
    }
    d1=\$(depth ${xls1})
    d2=\$(depth ${xls2})

    macs3 bdgdiff \\
        --t1 ${t1} --c1 ${c1} \\
        --t2 ${t2} --c2 ${c2} \\
        --d1 \$d1 --d2 \$d2 \\
        --o-prefix ${prefix} \\
        ${args}

    printf 'cond1\\t%s\\ncond2\\t%s\\nfactor\\t%s\\ndepth1_millions\\t%s\\ndepth2_millions\\t%s\\n' \\
        '${meta.cond1}' '${meta.cond2}' '${meta.factor}' "\$d1" "\$d2" > ${prefix}.contrast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version 2>&1 | sed 's/^macs3 //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cond1.bed ${prefix}_cond2.bed ${prefix}_common.bed
    printf 'cond1\\t%s\\ncond2\\t%s\\n' '${meta.cond1}' '${meta.cond2}' > ${prefix}.contrast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: stub
    END_VERSIONS
    """
}

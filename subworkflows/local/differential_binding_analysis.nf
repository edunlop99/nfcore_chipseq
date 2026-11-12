//
// Differential binding, one contrast at a time.
//
//   Group    = tissue + condition + treatment (e.g. PFC_WT_DMSO)
//   Contrast = two groups of one antibody that differ in exactly one factor.
//              The reference level (params.reference_*) is the denominator / cond1.
//
//   Both groups replicated (>= 2)  -> DiffBind
//   Both groups single-sample      -> MACS3 bdgdiff on the samples' own bedGraphs
//   One of each                    -> skipped with a warning
//

include { DIFFBIND      } from '../../../modules/local/diffbind.nf'
include { MACS3_BDGDIFF } from '../../../modules/nf-core/macs3/bdgdiff/main.nf'

workflow DIFFERENTIAL_BINDING {

    take:
    ch_peaks_bam_bai   // [ meta, peak, [ ip_bam, control_bam ], [ ip_bai, control_bai ] ]  control optional
    ch_bdg             // [ meta, [ *_treat_pileup.bdg, *_control_lambda.bdg ] ]  from MACS3_CALLPEAK, no --SPMR
    ch_xls             // [ meta, *_peaks.xls ]                                   from MACS3_CALLPEAK
    is_narrow_peak     // true = narrowPeak, false = broadPeak

    main:
    ch_versions = Channel.empty()

    def factors   = [ 'tissue', 'condition', 'treatment' ]
    def allowed   = params.contrast_factors ? params.contrast_factors.tokenize(',')*.trim() : factors
    def reference = [ tissue: params.reference_tissue, condition: params.reference_condition, treatment: params.reference_treatment ]

    //
    // All files for each sample, grouped per antibody.
    // sample = [ meta, peak, ip_bam, ip_bai, [ctrl_bam] or [], [ctrl_bai] or [], treat_bdg, ctrl_bdg, xls ]
    //
    ch_by_antibody = ch_peaks_bam_bai
        .map  { meta, peak, bams, bais -> [ meta.id, meta, peak, bams, bais ] }
        .join ( ch_bdg.map { meta, bdgs -> [ meta.id, bdgs ] }, failOnMismatch: true )
        .join ( ch_xls.map { meta, xls  -> [ meta.id, xls  ] }, failOnMismatch: true )
        .map  { id, meta, peak, bams, bais, bdgs, xls ->
            [ meta.antibody, [ meta, peak, bams[0], bais[0], bams.drop(1), bais.drop(1),
                               bdgs.find { it.name.endsWith('_treat_pileup.bdg') },
                               bdgs.find { it.name.endsWith('_control_lambda.bdg') },
                               xls ] ]
        }
        .groupTuple()

    //
    // Every contrast: [ contrast_meta, numerator_samples, denominator_samples ]
    //
    ch_contrasts = ch_by_antibody
        .flatMap { antibody, samples ->
            def by_group = samples.groupBy { s -> factors.collect { s[0][it] }.findAll { it }.join('_') }            
            def names    = by_group.keySet().sort()
            def out      = []
            for (int i = 0; i < names.size(); i++) {
                for (int j = i + 1; j < names.size(); j++) {
                    def a = by_group[names[i]][0][0]
                    def b = by_group[names[j]][0][0]
                    def differ = factors.findAll { a[it] != b[it] }
                    if (differ.size() != 1 || !(differ[0] in allowed)) continue

                    def f = differ[0]
                    def (den, num) = (b[f] == reference[f]) ? [ names[j], names[i] ] : [ names[i], names[j] ]
                    def dm = by_group[den][0][0]
                    def nm = by_group[num][0][0]
                    def fixed = factors.findAll { it != f }.collect { dm[it] }.findAll { it }.join('_') // if any factor is blank, leave out of the name
                    def cmeta = [
                        id          : "${antibody}_${f}_${nm[f]}_vs_${dm[f]}${fixed ? '_in_' + fixed : ''}".toString(),
                        antibody    : antibody,
                        factor      : f,
                        numerator   : num,
                        denominator : den,
                        cond1       : den,     // bdgdiff: *_cond1.bed = higher in reference
                        cond2       : num
                    ]
                    out << [ cmeta, by_group[num], by_group[den] ]
                }
            }
            if (!out) log.warn "[DIFFERENTIAL_BINDING] ${antibody}: no two groups differ in exactly one factor; no differential analysis"
            out
        }
        .branch { cmeta, num, den ->
            diffbind : num.size() >= 2 && den.size() >= 2
            bdgdiff  : num.size() == 1 && den.size() == 1
            skipped  : true
        }

    ch_contrasts.skipped.subscribe { cmeta, num, den ->
        log.warn "[DIFFERENTIAL_BINDING] Skipping ${cmeta.id}: ${num.size()} vs ${den.size()} replicates"
    }

    //
    // DiffBind: [ meta, metas, peaks, ip_bams, ip_bais, control_names, control_bams, control_bais ]
    //
    DIFFBIND (
        ch_contrasts.diffbind.map { cmeta, num, den ->
            def s = num + den
            [ cmeta,
              s.collect { it[0] }, s.collect { it[1] }, s.collect { it[2] }, s.collect { it[3] },
              s.collect { it[4] ? it[4][0].name : '' },
              s.collect { it[4] }.flatten().unique(),    // shared controls staged once
              s.collect { it[5] }.flatten().unique() ]
        },
        is_narrow_peak
    )
    ch_versions = ch_versions.mix(DIFFBIND.out.versions.first())

    //
    // bdgdiff: [ meta, t1, c1, xls1, t2, c2, xls2 ]  (1 = reference sample)
    //
    MACS3_BDGDIFF (
        ch_contrasts.bdgdiff.map { cmeta, num, den ->
            [ cmeta, den[0][6], den[0][7], den[0][8], num[0][6], num[0][7], num[0][8] ]
        }
    )
    ch_versions = ch_versions.mix(MACS3_BDGDIFF.out.versions.first())

    emit:
    diffbind_sheet   = DIFFBIND.out.samplesheet
    diffbind_results = DIFFBIND.out.results
    diffbind_rds     = DIFFBIND.out.rds
    bdgdiff_cond1    = MACS3_BDGDIFF.out.cond1    // higher in reference group
    bdgdiff_cond2    = MACS3_BDGDIFF.out.cond2    // higher in compared group
    bdgdiff_common   = MACS3_BDGDIFF.out.common
    versions         = ch_versions
}

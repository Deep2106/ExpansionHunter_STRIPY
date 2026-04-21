/*
========================================================================================
    SUBWORKFLOW: STRIPY_POSTPROCESS
========================================================================================
    Orchestrates all post-processing steps after ExpansionHunter:
      1. Fetch Stripy locus reference (once per pipeline run)
      2. Annotate per-sample VCFs locally (no VCF data sent externally)
      3. Build per-sample TSV reports    (/compare API — repeat counts only)
      4. Aggregate per-family Excel workbooks

    NOTE: REViewer visualisation is NOT included here.
    EH streaming mode (WES) does not produce realigned BAMs.
    Use the standalone utility:  scripts/run_reviewer.sh
----------------------------------------------------------------------------------------
*/

include { STRIPY_FETCH_LOCUS_REF        } from '../../modules/local/stripy/fetch_locus_ref'
include { STRIPY_ANNOTATE_VCF           } from '../../modules/local/stripy/annotate_vcf'
include { STRIPY_BUILD_SAMPLE_REPORT    } from '../../modules/local/stripy/build_sample_report'
include { STRIPY_AGGREGATE_COHORT_REPORT} from '../../modules/local/stripy/aggregate_cohort_report'

workflow STRIPY_POSTPROCESS {

    take:
    ch_vcf          // channel: [ val(meta), path(vcf) ]  from EXPANSIONHUNTER
    variant_catalog // path: variant catalog JSON
    reference       // path: reference FASTA
    ped             // path: PED file

    main:

    ch_versions = Channel.empty()

    // ── Step 1: Fetch Stripy locus reference (once) ───────────────────────
    STRIPY_FETCH_LOCUS_REF (
        variant_catalog
    )
    ch_versions = ch_versions.mix(STRIPY_FETCH_LOCUS_REF.out.versions)

    // ── Step 2: Annotate per-sample VCFs locally ─────────────────────────
    //    Uses local stripy_locus_reference.json — no VCF data sent externally
    STRIPY_ANNOTATE_VCF (
        ch_vcf,
        STRIPY_FETCH_LOCUS_REF.out.locus_ref
    )
    ch_versions = ch_versions.mix(STRIPY_ANNOTATE_VCF.out.versions.first())

    // ── Step 3: Collect all annotated VCFs for parent REPCN lookup ───────
    ch_all_annotated_vcfs = STRIPY_ANNOTATE_VCF.out.annotated_vcf
        .map { meta, vcf -> vcf }
        .collect()

    // ── Step 4: Build per-sample TSV reports ─────────────────────────────
    def hpo_file = file(params.hpo_file, checkIfExists: true)

    STRIPY_BUILD_SAMPLE_REPORT (
        STRIPY_ANNOTATE_VCF.out.annotated_vcf,
        STRIPY_FETCH_LOCUS_REF.out.locus_ref,
        hpo_file,
        ped,
        ch_all_annotated_vcfs
    )
    ch_versions = ch_versions.mix(STRIPY_BUILD_SAMPLE_REPORT.out.versions.first())

    // ── Step 5: Aggregate per-family Excel reports ────────────────────────
    ch_all_reports = STRIPY_BUILD_SAMPLE_REPORT.out.report
        .map { meta, tsv -> tsv }
        .collect()

    STRIPY_AGGREGATE_COHORT_REPORT (
        ch_all_reports,
        ped
    )
    ch_versions = ch_versions.mix(STRIPY_AGGREGATE_COHORT_REPORT.out.versions)

    emit:
    annotated_vcfs = STRIPY_ANNOTATE_VCF.out.annotated_vcf       // [ meta, vcf ]
    sample_reports = STRIPY_BUILD_SAMPLE_REPORT.out.report        // [ meta, tsv ]
    cohort_report  = STRIPY_AGGREGATE_COHORT_REPORT.out.cohort_excel
    versions       = ch_versions
}

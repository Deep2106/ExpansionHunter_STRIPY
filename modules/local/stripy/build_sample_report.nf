/*
========================================================================================
    STRIPY — BUILD SAMPLE REPORT
========================================================================================
    Generate per-sample STR annotation TSV by combining:
      - STRipy-annotated VCF
      - Stripy /locus reference JSON
      - HPO gene-to-phenotype file
      - PED file (family relationships)
      - All annotated VCFs directory (parent REPCN for DeNovo)
      - Stripy /compare API (called loci only)
    Output: per-sample TSV with all 48 required columns.
----------------------------------------------------------------------------------------
*/

process STRIPY_BUILD_SAMPLE_REPORT {
    tag        "${meta.id}"
    label      'process_low'

    container  params.containers.python

    publishDir "${params.outdir}/stripy/reports",
               mode: 'copy',
               pattern: "*.stripy_report.tsv"

    input:
    tuple val(meta), path(annotated_vcf)
    path locus_ref
    path hpo_file
    path ped
    path vcf_dir, stageAs: 'vcf_dir/*'   // all annotated VCFs staged into subdirectory
                                          // avoids filename collision with annotated_vcf input

    output:
    tuple val(meta), path("${meta.id}.stripy_report.tsv"), emit: report
    path "versions.yml",                                    emit: versions

    script:
    """
    python3 ${projectDir}/scripts/build_sample_report.py \\
        --vcf        ${annotated_vcf} \\
        --sample-id  ${meta.id} \\
        --locus-ref  ${locus_ref} \\
        --hpo-file   ${hpo_file} \\
        --ped        ${ped} \\
        --vcf-dir    vcf_dir \\
        --outdir     .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        stripy_compare_api: "https://api.stripy.org/compare"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.stripy_report.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        stripy_compare_api: "https://api.stripy.org/compare"
    END_VERSIONS
    """
}

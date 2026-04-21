/*
========================================================================================
    STRIPY — ANNOTATE VCF (LOCAL)
========================================================================================
    Annotates per-sample ExpansionHunter VCF using the pre-fetched local locus
    reference (stripy_locus_reference.json).

    Replaces the previous `curl POST https://api.stripy.org/annotateVCF` call.
    No patient VCF data is sent externally. Only /compare calls (repeat counts
    only, no patient identifiers) remain as external API interactions.

    Adds the same INFO fields as the API would:
        DISID    - disease identifier(s)
        DISNAME  - disease name(s)
        DISINHER - inheritance pattern(s)
        DISRANGE - normal range Min|Max per disease
----------------------------------------------------------------------------------------
*/

process STRIPY_ANNOTATE_VCF {
    tag        "${meta.id}"
    label      'process_low'

    container  params.containers.python

    publishDir "${params.outdir}/stripy/annotated_vcf",
               mode: 'copy',
               pattern: "*.Reannotated.vcf"

    input:
    tuple val(meta), path(vcf)
    path locus_ref

    output:
    tuple val(meta), path("${meta.id}.Reannotated.vcf"), emit: annotated_vcf
    path "versions.yml",                                  emit: versions

    script:
    """
    python3 ${projectDir}/scripts/annotate_vcf_local.py \\
        --vcf       ${vcf} \\
        --locus-ref ${locus_ref} \\
        --output    ${meta.id}.Reannotated.vcf

    # Verify output is a valid VCF
    if ! grep -q "^#CHROM" ${meta.id}.Reannotated.vcf; then
        echo "[ERROR] Local annotation produced invalid VCF for ${meta.id}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        annotate_vcf_local: "local — no external API call"
    END_VERSIONS
    """

    stub:
    """
    cp ${vcf} ${meta.id}.Reannotated.vcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11.0"
        annotate_vcf_local: "local — no external API call"
    END_VERSIONS
    """
}

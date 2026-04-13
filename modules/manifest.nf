/*
========================================================================================
    MANIFEST MODULE
========================================================================================
    Write a TSV manifest of all completed samples after ExpansionHunter finishes.
    One row per sample: sample_id, family_id, sex, phenotype, json_path.
    Used for tracking pipeline completion and downstream QC.
----------------------------------------------------------------------------------------
*/

process WRITE_MANIFEST {
    tag   "manifest"
    label 'process_single'

    container params.containers.python

    publishDir "${params.outdir}/pipeline_info",
               mode: params.publish_dir_mode

    input:
    val rows   // collected list of TSV row strings from EXPANSIONHUNTER.out.json

    output:
    path "completed_samples.tsv", emit: manifest
    path "versions.yml",          emit: versions

    script:
    """
    python3 - <<'PYEOF'
import sys
from pathlib import Path

rows = """${rows.join('\n')}""".strip().splitlines()

out = Path("completed_samples.tsv")
with open(out, 'w') as f:
    f.write("sample_id\\tfamily_id\\tsex\\tphenotype\\tjson_path\\n")
    for row in rows:
        if row.strip():
            f.write(row.strip() + "\\n")

print(f"[INFO] Manifest written : {out}  ({len(rows)} samples)")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch completed_samples.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.14.2
    END_VERSIONS
    """
}

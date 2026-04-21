#!/usr/bin/env python3
"""
annotate_vcf_local.py
=====================
Local replacement for the STRipy POST /annotateVCF API call.

Reads the pre-fetched stripy_locus_reference.json and annotates an
ExpansionHunter VCF with the same INFO fields that /annotateVCF would add:

    DISID    - comma-separated disease IDs
    DISNAME  - comma-separated disease names (spaces replaced with underscores)
    DISINHER - comma-separated inheritance patterns
    DISRANGE - comma-separated NormalMin|NormalMax per disease

Loci not found in the reference are passed through unchanged (no DISID added),
which exactly mirrors the /annotateVCF behaviour.

Usage:
    python3 annotate_vcf_local.py \\
        --vcf         sample.vcf \\
        --locus-ref   stripy_locus_reference.json \\
        --output      sample.Reannotated.vcf
"""

import argparse
import json
import re
import sys
from pathlib import Path


# ── New VCF INFO header lines to inject ──────────────────────────────────────

NEW_INFO_HEADERS = [
    '##INFO=<ID=DISID,Number=.,Type=String,Description="STRipy disease identifier(s)">',
    '##INFO=<ID=DISNAME,Number=.,Type=String,Description="STRipy disease name(s); spaces replaced with underscores">',
    '##INFO=<ID=DISINHER,Number=.,Type=String,Description="STRipy disease inheritance pattern(s)">',
    '##INFO=<ID=DISRANGE,Number=.,Type=String,Description="STRipy normal range as Min|Max per disease">',
]


# ── Locus reference helpers ───────────────────────────────────────────────────

def load_locus_reference(ref_path: str) -> dict:
    """
    Load stripy_locus_reference.json and index by normalised locus ID.
    Returns dict keyed by locus ID (e.g. 'HTT', 'RFC1', 'HOXA13_1').
    """
    p = Path(ref_path)
    if not p.exists():
        sys.exit(f"[ERROR] Locus reference not found: {ref_path}")
    with open(p) as f:
        data = json.load(f)
    ref = {}
    for entry in data:
        if isinstance(entry, dict) and 'Locus' in entry:
            ref[entry['Locus']] = entry
    print(f"[INFO] Loaded {len(ref)} loci from reference", file=sys.stderr)
    return ref


def normalise_locus_id(repid: str) -> str:
    """
    Normalise REPID to a plain locus ID for reference lookup.

    EH VCF REPID field can contain motif suffixes e.g. 'RFC1:AAGGG'.
    The locus reference uses the bare locus ID ('RFC1').
    """
    return repid.split(':')[0]


def build_disease_annotations(diseases: dict) -> tuple[str, str, str, str]:
    """
    Build comma-separated DISID, DISNAME, DISINHER, DISRANGE strings
    from the Diseases dict of a locus reference entry.

    Returns (dis_id, dis_name, dis_inher, dis_range) — all comma-joined.
    """
    dis_ids    = []
    dis_names  = []
    dis_inher  = []
    dis_ranges = []

    for dis_id, dis in diseases.items():
        dis_ids.append(dis_id)

        name = dis.get('DiseaseName', dis_id).replace(' ', '_')
        dis_names.append(name)

        inher = dis.get('Inheritance', 'Unknown')
        # Normalise to short form to match /annotateVCF output
        # e.g. "Autosomal dominant" → "AD", "Autosomal recessive" → "AR"
        inher_map = {
            'Autosomal dominant':  'AD',
            'Autosomal recessive': 'AR',
            'X-linked dominant':   'XD',
            'X-linked recessive':  'XR',
            'X-linked':            'XL',
            'Mitochondrial':       'MT',
        }
        dis_inher.append(inher_map.get(inher, inher))

        nr = dis.get('NormalRange') or {}
        nr_min = nr.get('Min', 'Unknown')
        nr_max = nr.get('Max', 'Unknown')
        dis_ranges.append(f"{nr_min}|{nr_max}")

    return (
        ','.join(dis_ids),
        ','.join(dis_names),
        ','.join(dis_inher),
        ','.join(dis_ranges),
    )


# ── VCF INFO field helpers ────────────────────────────────────────────────────

def parse_info(info_str: str) -> dict:
    """Parse semicolon-delimited INFO string into dict."""
    info = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            info[k] = v
        else:
            info[field] = True
    return info


def serialize_info(info: dict) -> str:
    """Serialize INFO dict back to semicolon-delimited string."""
    parts = []
    for k, v in info.items():
        if v is True:
            parts.append(k)
        else:
            parts.append(f"{k}={v}")
    return ';'.join(parts)


# ── Main annotation logic ─────────────────────────────────────────────────────

def annotate_vcf(vcf_path: str, locus_ref: dict, output_path: str) -> None:
    n_annotated = 0
    n_skipped   = 0
    n_records   = 0
    header_written = False

    with open(vcf_path) as fh_in, open(output_path, 'w') as fh_out:
        for line in fh_in:
            # ── Header lines ─────────────────────────────────────────────
            if line.startswith('#'):
                # Inject new INFO headers just before the #CHROM line
                if line.startswith('#CHROM') and not header_written:
                    for hdr in NEW_INFO_HEADERS:
                        fh_out.write(hdr + '\n')
                    header_written = True
                fh_out.write(line)
                continue

            # ── Data lines ───────────────────────────────────────────────
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                fh_out.write(line)
                continue

            n_records += 1
            info = parse_info(cols[7])

            # Get REPID — required to identify the locus
            repid = info.get('REPID')
            if not repid:
                fh_out.write(line)
                n_skipped += 1
                continue

            locus_id = normalise_locus_id(repid)
            entry    = locus_ref.get(locus_id)

            if not entry or not entry.get('Diseases'):
                # Locus not in reference or no disease data — pass through unchanged
                fh_out.write(line)
                n_skipped += 1
                continue

            # Build disease annotation fields
            dis_id, dis_name, dis_inher, dis_range = build_disease_annotations(
                entry['Diseases']
            )

            # Inject into INFO (append after existing fields, same as /annotateVCF)
            info['DISID']    = dis_id
            info['DISNAME']  = dis_name
            info['DISINHER'] = dis_inher
            info['DISRANGE'] = dis_range

            cols[7] = serialize_info(info)
            fh_out.write('\t'.join(cols) + '\n')
            n_annotated += 1

    print(
        f"[INFO] Records: {n_records} total | "
        f"{n_annotated} annotated | {n_skipped} skipped (not in reference)",
        file=sys.stderr
    )


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Local STRipy VCF annotation (replaces POST /annotateVCF)'
    )
    p.add_argument('--vcf',       required=True, help='Input ExpansionHunter VCF')
    p.add_argument('--locus-ref', required=True, help='stripy_locus_reference.json')
    p.add_argument('--output',    required=True, help='Output annotated VCF path')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    locus_ref = load_locus_reference(args.locus_ref)
    annotate_vcf(args.vcf, locus_ref, args.output)
    print(f"[INFO] Annotated VCF written to: {args.output}", file=sys.stderr)

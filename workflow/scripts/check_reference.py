"""Compare CRAM ``@SQ`` contigs to a SAM sequence dictionary.

CRAM decoding needs the same reference the aligner used. Broad/AnVIL CRAMs
are expected to match ``Homo_sapiens_assembly38`` (chr-prefixed GRCh38), but
the FASTA is a parameter: this check fails when MD5s disagree.
"""

from __future__ import annotations

import argparse
import json
import sys


def contigs_from_sam_header(text: str) -> dict[str, str | None]:
    """Map contig name to MD5 (``M5``) from ``@SQ`` lines."""
    found: dict[str, str | None] = {}
    for line in text.splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = {}
        for part in line.split("\t")[1:]:
            if ":" in part:
                key, value = part.split(":", 1)
                fields[key] = value
        name = fields.get("SN")
        if name:
            found[name] = fields.get("M5")
    return found


def compare_headers(cram_header: str, reference_dict: str) -> dict:
    """Return a JSON-serialisable report. ``ok`` is false when the CRAM cannot be decoded."""
    cram = contigs_from_sam_header(cram_header)
    reference = contigs_from_sam_header(reference_dict)
    missing = []
    md5_mismatch = []
    md5_missing = []
    for name, md5 in cram.items():
        if name not in reference:
            missing.append(name)
            continue
        ref_md5 = reference[name]
        if not md5 or not ref_md5:
            md5_missing.append(name)
        elif md5.lower() != ref_md5.lower():
            md5_mismatch.append({"contig": name, "cram_md5": md5, "reference_md5": ref_md5})
    prefix_hint = None
    if missing and reference:
        cram_names = set(cram)
        ref_names = set(reference)
        if {f"chr{name}" for name in cram_names} <= ref_names or {name[3:] for name in cram_names if name.startswith("chr")} <= ref_names:
            prefix_hint = "Contig names differ by a chr prefix. Broad assembly38 uses chr1..chrM."
    ok = not missing and not md5_mismatch
    return {
        "ok": ok,
        "cram_contigs": len(cram),
        "reference_contigs": len(reference),
        "missing_from_reference": missing,
        "md5_mismatch": md5_mismatch,
        "md5_not_in_both_headers": md5_missing,
        "prefix_hint": prefix_hint,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cram-header", required=True, help="samtools view -H output")
    parser.add_argument("--reference-dict", required=True, help="samtools dict / Picard dict output")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    report = compare_headers(
        open(args.cram_header, encoding="utf-8", errors="replace").read(),
        open(args.reference_dict, encoding="utf-8", errors="replace").read(),
    )
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    if not report["ok"]:
        print(json.dumps(report), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

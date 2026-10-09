"""Memory-safe VCF access.

The joint SNV VCF is about 20 GB. Callers pass a BED and this module fetches
those regions with bcftools, tabix, pysam, or cyvcf2 when an index exists.
Files larger than ``STREAM_MAX_BYTES`` are never read end to end.
"""

from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

STREAM_MAX_BYTES = 512 * 1024 * 1024

_CONTIG_RE = re.compile(r"ID=([^,>]+)")
_INTERVAL_RE = re.compile(r"(chr)?([0-9XYM]+|MT):(\d+)-(\d+)", re.I)


def _open_text(path: Path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def chrom_bare(chrom: str) -> str:
    text = str(chrom)
    if text.lower().startswith("chr"):
        text = text[3:]
    return "M" if text.upper() in {"MT", "M"} else text


def index_path(path: Path) -> Path | None:
    for suffix in (".tbi", ".csi"):
        candidate = Path(str(path) + suffix)
        if candidate.exists():
            return candidate
    return None


def read_samples(path: Path) -> list[str]:
    """Sample IDs from the ``#CHROM`` line. Stops at the header, so a 20 GB VCF stays on disk."""
    with _open_text(Path(path)) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"No #CHROM header in {path}")


def read_contigs(path: Path) -> set[str]:
    contigs: set[str] = set()
    with _open_text(Path(path)) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                break
            if line.startswith("##contig"):
                match = _CONTIG_RE.search(line)
                if match:
                    contigs.add(match.group(1))
    return contigs


def style_chrom(chrom: str, contigs: set[str]) -> str:
    """Return ``chrom`` spelled the way ``contigs`` spells it."""
    if not contigs:
        return str(chrom)
    bare = chrom_bare(chrom)
    with_chr = "chrM" if bare == "M" else f"chr{bare}"
    if with_chr in contigs:
        return with_chr
    if bare in contigs:
        return bare
    return str(chrom)


def parse_info(info: str) -> dict[str, str | bool]:
    out: dict[str, str | bool] = {}
    if not info or info == ".":
        return out
    for part in info.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
        else:
            out[part] = True
    return out


def parse_format(fmt: str, sample: str) -> dict[str, str]:
    if not fmt or fmt == ".":
        return {}
    keys = fmt.split(":")
    vals = sample.split(":") if sample else []
    return {key: (vals[i] if i < len(vals) else ".") for i, key in enumerate(keys)}


def parse_record(line: str, samples: list[str]) -> dict:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 8:
        raise ValueError(f"VCF record has {len(parts)} columns")
    chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
    fmt = parts[8] if len(parts) > 8 else ""
    sample_fields = parts[9:]
    parsed = {
        name: parse_format(fmt, sample_fields[i] if i < len(sample_fields) else ".")
        for i, name in enumerate(samples)
    }
    return {
        "chrom": chrom,
        "pos": int(pos),
        "id": vid,
        "ref": ref,
        "alt": alt,
        "qual": qual,
        "filter": filt,
        "info": parse_info(info),
        "format": fmt,
        "samples": parsed,
    }


class IntervalIndex:
    """Half-open BED intervals grouped by chromosome (``chr`` prefix stripped)."""

    def __init__(self, records: list[tuple[str, int, int, dict]]):
        self.by_chrom: dict[str, list[tuple[int, int, dict]]] = {}
        original: list[tuple[str, int, int]] = []
        for chrom, start, end, payload in records:
            self.by_chrom.setdefault(chrom_bare(chrom), []).append((int(start), int(end), payload))
            original.append((str(chrom), int(start), int(end)))
        for chrom in self.by_chrom:
            self.by_chrom[chrom].sort()
        self.merged = _merge_regions(original)

    @classmethod
    def from_bed(cls, path: str | Path, flank: int = 0) -> "IntervalIndex":
        records: list[tuple[str, int, int, dict]] = []
        with Path(path).open(encoding="utf-8") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                start = max(int(parts[1]) - flank, 0)
                end = int(parts[2]) + flank
                payload = {
                    "chrom": parts[0],
                    "start": start,
                    "end": end,
                    "gene_id": parts[3] if len(parts) > 3 else "",
                    "gene": (parts[4] if len(parts) > 4 else "").upper(),
                }
                records.append((parts[0], start, end, payload))
        return cls(records)

    def overlap(self, chrom: str, start: int, end: int) -> list[dict]:
        hits: list[dict] = []
        for iv_start, iv_end, payload in self.by_chrom.get(chrom_bare(chrom), []):
            if iv_start >= end:
                break
            if iv_end > start:
                hits.append(payload)
        return hits

    def write_bed(self, path: str | Path, contigs: set[str] | None = None) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        contigs = contigs or set()
        with path.open("w", encoding="utf-8", newline="\n") as fh:
            for chrom, start, end in self.merged:
                fh.write(f"{style_chrom(chrom, contigs)}\t{start}\t{end}\n")
        return path


def _merge_regions(regions: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    grouped: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in regions:
        grouped.setdefault(chrom, []).append((start, end))
    merged: list[tuple[str, int, int]] = []
    for chrom, spans in grouped.items():
        spans.sort()
        cur_s, cur_e = spans[0]
        for start, end in spans[1:]:
            if start <= cur_e:
                cur_e = max(cur_e, end)
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = start, end
        merged.append((chrom, cur_s, cur_e))
    return merged


def ensure_tabix(path: Path) -> bool:
    """Index a bgzip VCF when ``tabix`` or ``bcftools`` is on ``PATH``. Return whether an index exists."""
    path = Path(path)
    if index_path(path):
        return True
    if path.stat().st_size > STREAM_MAX_BYTES:
        return False
    if not str(path).endswith((".gz", ".bgz")):
        return False
    tabix = shutil.which("tabix")
    bcftools = shutil.which("bcftools")
    try:
        if tabix:
            subprocess.run([tabix, "-p", "vcf", str(path)], check=True, capture_output=True, text=True)
        elif bcftools:
            subprocess.run([bcftools, "index", "-t", str(path)], check=True, capture_output=True, text=True)
        else:
            return False
    except subprocess.CalledProcessError:
        return False
    return index_path(path) is not None


def _which_backend(path: Path) -> str | None:
    if index_path(path) is None:
        return None
    if shutil.which("bcftools"):
        return "bcftools"
    if shutil.which("tabix"):
        return "tabix"
    try:
        import pysam  # noqa: F401

        return "pysam"
    except ImportError:
        pass
    try:
        import cyvcf2  # noqa: F401

        return "cyvcf2"
    except ImportError:
        return None


def iter_vcf(
    path: str | Path,
    regions: str | Path | None = None,
    *,
    allow_full_scan: bool | None = None,
):
    """Yield parsed records.

    ``regions`` is a BED path. Without it, a file larger than ``STREAM_MAX_BYTES``
    raises ``RuntimeError`` instead of being scanned.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    size = path.stat().st_size
    if regions is None:
        allowed = allow_full_scan if allow_full_scan is not None else size <= STREAM_MAX_BYTES
        if not allowed:
            raise RuntimeError(
                f"Refusing to scan {path} ({size} bytes) without a region BED. "
                "Query the panel BED with bcftools or tabix."
            )
        yield from _iter_python(path, None)
        return

    backend = _which_backend(path)
    if backend:
        yield from _iter_backend(path, Path(regions), backend)
        return
    if size > STREAM_MAX_BYTES and allow_full_scan is not True:
        raise RuntimeError(
            f"{path} is larger than {STREAM_MAX_BYTES} bytes and has no tabix index. "
            "Install bcftools or tabix and index the VCF. The joint SNV file must be queried by region."
        )
    yield from _iter_python(path, IntervalIndex.from_bed(regions))


def _iter_python(path: Path, index: IntervalIndex | None):
    samples: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if not line.strip():
                continue
            record = parse_record(line, samples)
            if index is not None and not _record_overlaps(record, index):
                continue
            yield record


def _record_overlaps(record: dict, index: IntervalIndex) -> bool:
    start0 = record["pos"] - 1
    end_info = record["info"].get("END")
    try:
        end0 = int(end_info) if isinstance(end_info, str) else record["pos"]
    except ValueError:
        end0 = record["pos"]
    if end0 < start0:
        end0 = record["pos"]
    if index.overlap(record["chrom"], start0, end0):
        return True
    raw = record["info"].get("CPX_INTERVALS")
    if isinstance(raw, str):
        for match in _INTERVAL_RE.finditer(raw):
            chrom = ("chr" if match.group(1) else "") + match.group(2)
            if index.overlap(chrom, int(match.group(3)) - 1, int(match.group(4))):
                return True
    return False


def _iter_backend(path: Path, regions: Path, backend: str):
    samples = read_samples(path)
    contigs = read_contigs(path)
    index = IntervalIndex.from_bed(regions)
    with tempfile.TemporaryDirectory() as tmp:
        query = index.write_bed(Path(tmp) / "regions.bed", contigs)
        if backend == "bcftools":
            lines = _stream_command(["bcftools", "view", "-R", str(query), "-H", str(path)])
        elif backend == "tabix":
            lines = _stream_command(["tabix", "-R", str(query), str(path)])
        elif backend == "pysam":
            lines = _pysam_lines(path, index, contigs)
        else:
            lines = _cyvcf2_lines(path, index, contigs)
        for line in lines:
            if not line or line.startswith("#"):
                continue
            yield parse_record(line, samples)


def _stream_command(cmd: list[str]):
    """Yield stdout lines. The VCF itself is not read into memory."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            yield line.rstrip("\n")
    finally:
        proc.stdout.close()
        code = proc.wait()
    if code != 0:
        raise RuntimeError(f"{cmd[0]} exited {code} while querying {cmd[-1]}")


def _pysam_lines(path: Path, index: IntervalIndex, contigs: set[str]):
    import pysam

    tabix = pysam.TabixFile(str(path))
    try:
        for chrom, start, end in index.merged:
            named = style_chrom(chrom, contigs)
            try:
                fetched = tabix.fetch(named, start, end)
            except ValueError:
                continue
            yield from fetched
    finally:
        tabix.close()


def _cyvcf2_lines(path: Path, index: IntervalIndex, contigs: set[str]):
    from cyvcf2 import VCF

    vcf = VCF(str(path))
    for chrom, start, end in index.merged:
        named = style_chrom(chrom, contigs)
        # cyvcf2 regions are 1-based and inclusive.
        for rec in vcf(f"{named}:{start + 1}-{end}"):
            yield str(rec)

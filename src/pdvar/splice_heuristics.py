"""Sullivan et al. 2025 donor and acceptor splice-disruption heuristics.

Figures 3 and 4 of Sullivan et al., Am J Hum Genet 2025;112:764-778
(doi:10.1016/j.ajhg.2025.02.012). A matching SNV is labelled with the
published subgroup, %SAV, and n. Other alleles are left blank. The labels
do not change the evidence score.

Positions are transcript-oriented: ``+`` is the donor intron, ``-`` is the
acceptor intron, ``E-`` is the exon end, and ``E+1`` is the first exonic base.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .normalize import apply_gene_aliases
from .vcf_stream import chrom_bare

COLUMNS = ("splice_rule", "splice_subgroup", "splice_sav_pct", "splice_heuristic_n")

# Published spliceogenicity and n (CCICB/SplicingHeuristics2025 heuristic_values.tsv).
_PUBLISHED: dict[tuple[str, str], tuple[float, int]] = {
    ("DD1", "context"): (81.2, 16),
    ("DD1", "standard"): (100.0, 28),
    ("DD2", "standard"): (99.9, 2218),
    ("DD3", "standard"): (100.0, 78),
    ("DD4", "standard"): (99.7, 583),
    ("DD5", "standard"): (95.4, 108),
    ("DD6", "auxiliary"): (28.6, 14),
    ("DD6", "standard"): (95.1, 446),
    ("DD7", "auxiliary"): (25.0, 12),
    ("DD7", "standard"): (94.0, 667),
    ("DD8", "auxiliary"): (25.0, 20),
    ("DD8", "standard"): (93.6, 344),
    ("DD9", "auxiliary"): (19.0, 63),
    ("DD9", "standard"): (91.7, 120),
    ("DD10", "auxiliary"): (14.3, 28),
    ("DD10", "standard"): (91.5, 165),
    ("DD11", "auxiliary"): (38.6, 57),
    ("DD11", "standard"): (91.1, 146),
    ("DD12", "auxiliary"): (45.3, 64),
    ("DD12", "standard"): (68.2, 22),
    ("DA1", "auxiliary"): (28.0, 232),
    ("DA1", "context"): (78.0, 100),
    ("DA1", "standard"): (97.7, 304),
    ("DA2", "standard"): (99.7, 2817),
    ("DA3", "context"): (70.4, 81),
    ("DA3", "standard"): (99.2, 124),
    ("DA4", "auxiliary"): (17.4, 23),
    ("DA4", "context"): (23.8, 21),
    ("DA4", "standard"): (96.8, 157),
    ("DA5", "auxiliary"): (43.3, 164),
    ("DA5", "standard"): (93.1, 102),
    ("DA6", "auxiliary"): (14.8, 27),
    ("DA6", "standard"): (66.7, 96),
    ("DA7", "auxiliary"): (3.9, 950),
    ("DA7", "context"): (19.2, 151),
    ("DA7", "standard"): (38.9, 18),
    ("DA8", "context"): (7.7, 39),
    ("DA8", "standard"): (33.3, 57),
    ("DA9", "auxiliary"): (4.6, 524),
    ("DA9", "context"): (3.3, 153),
    ("DA9", "standard"): (25.4, 508),
}

# Most frequent base at E-3..+6. Figure 1 writes the donor consensus as AG|GTRAG;
# +3 accepts either purine. The GC row uses the compensatory GC-donor motif
# (strong G at E-1 and +5, C at +2).
_GT_IDEAL = {-3: "C", -2: "A", -1: "G", 1: "G", 2: "T", 3: "AG", 4: "A", 5: "G", 6: "T"}
_GC_IDEAL = {-3: "C", -2: "A", -1: "G", 1: "G", 2: "C", 3: "A", 4: "A", 5: "G", 6: "T"}
_DONOR_OFFSETS = (-3, -2, -1, 1, 2, 3, 4, 5, 6)
_ACCEPTOR_OFFSETS = tuple(range(-52, 0)) + (1,)
_HGVS_OFFSET = re.compile(r"c\.(?:\d+)(?P<sign>[+-])(?P<dist>\d+)")
_HGVS_TRANSCRIPT = re.compile(r"^(?P<id>(?:ENST|NM_|NR_|XM_|XR_)\d+)")
_GTF_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')
_PYRIMIDINE = frozenset("CT")


@dataclass(frozen=True)
class Exon:
    chrom: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive
    strand: str
    gene: str
    transcript: str


class MemoryFasta:
    """1-based inclusive fetch over an in-memory reference, for tests."""

    def __init__(self, sequences: dict[str, str]):
        self.sequences = {chrom_bare(name): value.upper() for name, value in sequences.items()}

    def fetch(self, chrom: str, start: int, end: int) -> str:
        seq = self.sequences[chrom_bare(chrom)]
        if start < 1 or end > len(seq) or start > end:
            raise ValueError(f"interval {chrom}:{start}-{end} is outside the reference")
        return seq[start - 1 : end]


class IndexedFasta:
    """Random access via a ``.fai`` index. Line wrapping is honoured."""

    def __init__(self, fasta: str | Path, index: str | Path):
        self._path = Path(fasta)
        self._entries: dict[str, tuple[int, int, int, int]] = {}
        for line in Path(index).read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            name, length, offset, line_bases, line_width = line.split("\t")[:5]
            self._entries[chrom_bare(name)] = (int(length), int(offset), int(line_bases), int(line_width))
        self._handle = self._path.open("rb")

    def fetch(self, chrom: str, start: int, end: int) -> str:
        length, offset, line_bases, line_width = self._entries[chrom_bare(chrom)]
        if start < 1 or end > length or start > end:
            raise ValueError(f"interval {chrom}:{start}-{end} is outside the reference")
        parts: list[str] = []
        pos = start
        while pos <= end:
            zero = pos - 1
            byte = offset + (zero // line_bases) * line_width + (zero % line_bases)
            take = min(end - pos + 1, line_bases - (zero % line_bases))
            self._handle.seek(byte)
            parts.append(self._handle.read(take).decode("ascii"))
            pos += take
        return "".join(parts).upper()

    def close(self) -> None:
        self._handle.close()


def load_gtf_exons(
    path: str | Path,
    genes: set[str] | None = None,
    aliases: dict[str, str] | None = None,
) -> list[Exon]:
    """Exon rows with strand. ``genes`` is matched after aliasing and upper-casing."""
    wanted = {gene.upper() for gene in genes} if genes is not None else None
    alias_map = {key.upper(): value.upper() for key, value in (aliases or {}).items()}
    opener = gzip.open if str(path).endswith(".gz") else open
    exons: list[Exon] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = dict(_GTF_ATTR.findall(parts[8]))
            gene = (attrs.get("gene_name") or attrs.get("gene_id") or "").upper()
            gene = alias_map.get(gene, gene)
            if wanted is not None and gene not in wanted:
                continue
            strand = parts[6]
            if strand not in "+-":
                continue
            transcript = (attrs.get("transcript_id") or "").split(".")[0]
            exons.append(
                Exon(
                    chrom=chrom_bare(parts[0]),
                    start=int(parts[3]) - 1,
                    end=int(parts[4]),
                    strand=strand,
                    gene=gene,
                    transcript=transcript,
                )
            )
    return exons


def _blank(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["splice_rule"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df["splice_subgroup"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df["splice_sav_pct"] = pd.Series(pd.NA, index=df.index, dtype="Float64")
    df["splice_heuristic_n"] = pd.Series(pd.NA, index=df.index, dtype="Int64")
    return df


def _more_preferred(ref: str, alt: str, order: str) -> bool:
    """True when ``ref`` outranks ``alt``. ``=`` ties a pair (``G>A=C>T``)."""
    rank: dict[str, int] = {}
    score = len(order)
    rank[order[0]] = score
    index = 1
    while index < len(order):
        operator = order[index]
        base = order[index + 1]
        if operator == ">":
            score -= 1
        rank[base] = score
        index += 2
    return rank[ref] > rank[alt]


def _with_alt(bases: dict[int, str], offset: int, alt: str) -> dict[int, str]:
    updated = dict(bases)
    updated[offset] = alt
    return updated


def _ideal_count(bases: dict[int, str], ideal: dict[int, str]) -> int:
    return sum(bases[offset] in ideal[offset] for offset in _DONOR_OFFSETS)


def _dinuc(bases: dict[int, str], left: int, right: int) -> str:
    return bases[left] + bases[right]


def _ppt_strong(bases: dict[int, str]) -> bool:
    """Twelve consecutive pyrimidines, or 17 pyrimidines, from -24 to -4."""
    run = 0
    longest = 0
    total = 0
    for offset in range(-24, -3):
        if bases.get(offset) in _PYRIMIDINE:
            total += 1
            run += 1
            longest = max(longest, run)
        else:
            run = 0
    return longest >= 12 or total >= 17


def _bps_positions(bases: dict[int, str]) -> list[int]:
    """Branchpoint A's (YNA) whose A sits between -50 and -17."""
    found = []
    for aden in range(-50, -16):
        if bases.get(aden) == "A" and bases.get(aden - 2) in _PYRIMIDINE:
            found.append(aden)
    return found


def _created_ag_as(ref: dict[int, str], alt: dict[int, str], changed: int) -> list[int]:
    """Positions of the A in an AG that the variant creates."""
    created = []
    for aden in (changed, changed - 1):
        if alt.get(aden) == "A" and alt.get(aden + 1) == "G":
            if not (ref.get(aden) == "A" and ref.get(aden + 1) == "G"):
                created.append(aden)
    return created


def _call(rule: str, subgroup: str) -> tuple[str, str]:
    return rule, subgroup


def classify_donor(key: str, ref: dict[int, str], alt_base: str) -> tuple[str, str] | None:
    """One Figure 3 terminal, or None when the position is outside the tree."""
    try:
        changed = -int(key[2:]) if key.startswith("E-") else int(key)
        if any(ref[offset] not in "ACGT" for offset in _DONOR_OFFSETS):
            return None
    except (KeyError, ValueError):
        return None
    alt = _with_alt(ref, changed, alt_base)

    if changed in (1, 2):
        dinuc = _dinuc(ref, 1, 2)
        if dinuc == "GT":
            return _call("DD2", "standard")
        if dinuc == "GC":
            if _dinuc(alt, 1, 2) == "GT" or _ideal_count(alt, _GC_IDEAL) >= 8:
                return _call("DD1", "context")
            return _call("DD1", "standard")
        return None

    if changed == -3:
        if alt_base != "T":
            return _call("DD12", "auxiliary")
        if ref[-2] == "A":
            return _call("DD3", "standard")
        return _call("DD12", "standard")

    if changed == -2:
        if ref[-3] == "T" and alt_base == "A":
            return _call("DD3", "standard")
        if ref[-2] == "A":
            return _call("DD11", "standard")
        return _call("DD11", "auxiliary")

    if changed == -1:
        if ref[-1] == "G" and ref[5] != "G":
            return _call("DD4", "standard")
        if _more_preferred(ref[-1], alt_base, "G>A>C>T"):
            return _call("DD7", "standard")
        return _call("DD7", "auxiliary")

    if changed == 5:
        if ref[5] != "G":
            return _call("DD6", "auxiliary")
        if ref[-1] == "G":
            return _call("DD6", "standard")
        return _call("DD4", "standard")

    if changed in (3, 4) and _dinuc(alt, 3, 4) in {"CC", "GT"}:
        return _call("DD5", "standard")

    if changed == 3:
        if _more_preferred(ref[3], alt_base, "A>G>C=T"):
            return _call("DD8", "standard")
        return _call("DD8", "auxiliary")

    if changed == 4:
        if ref[4] == "A":
            return _call("DD10", "standard")
        return _call("DD10", "auxiliary")

    if changed == 6:
        if ref[6] == "T":
            return _call("DD9", "standard")
        return _call("DD9", "auxiliary")
    return None


def classify_acceptor(key: str, ref: dict[int, str], alt_base: str) -> tuple[str, str] | None:
    """One Figure 4 terminal, or None when the position is outside the tree."""
    try:
        changed = 1 if key == "E+1" else int(key)
    except ValueError:
        return None
    if changed == 0 or ref.get(changed) not in "ACGT":
        return None
    alt = _with_alt(ref, changed, alt_base)

    if changed in (-1, -2) and ref.get(-2) == "A" and ref.get(-1) == "G":
        return _call("DA2", "standard")

    created = _created_ag_as(ref, alt, changed)
    if created:
        if any(-13 <= aden <= -6 for aden in created):
            return _call("DA1", "standard")
        branchpoints = _bps_positions(ref)
        for aden in created:
            if aden > -14:
                continue
            upstream = [point for point in branchpoints if point < aden]
            if upstream and aden - max(upstream) >= 3:
                return _call("DA1", "context")
        return _call("DA1", "auxiliary")

    if changed == 1:
        if ref.get(1) == "G" and ref.get(-3) != "C":
            if ref.get(-3) == "T":
                return _call("DA3", "context")
            return _call("DA3", "standard")
        if ref.get(1) in "ACGT" and _more_preferred(ref[1], alt_base, "G>A=C>T"):
            return _call("DA6", "standard")
        return _call("DA6", "auxiliary")

    if changed == -3:
        if ref.get(-3) == "C" and ref.get(1) == "G" and alt_base == "T":
            return _call("DA4", "context")
        if ref.get(-3) == "C" and ref.get(1) == "G":
            return _call("DA4", "standard")
        if ref.get(-3) in "ACGT" and _more_preferred(ref[-3], alt_base, "C>T>A>G"):
            return _call("DA4", "standard")
        return _call("DA4", "auxiliary")

    if changed == -4:
        if ref.get(-6) == "T" and ref.get(-5) == "T" or _ppt_strong(ref):
            return _call("DA8", "context")
        return _call("DA8", "standard")

    if changed in (-5, -6):
        if ref.get(changed) == "T" and alt_base != "C":
            return _call("DA5", "standard")
        return _call("DA5", "auxiliary")

    if -24 <= changed <= -7:
        if _ppt_strong(ref):
            return _call("DA9", "context")
        if alt_base == "T" or ref.get(changed) not in _PYRIMIDINE:
            return _call("DA9", "auxiliary")
        return _call("DA9", "standard")

    if -50 <= changed <= -25:
        ref_n = len(_bps_positions(ref))
        alt_n = len(_bps_positions(alt))
        if alt_n < ref_n and alt_n == 0:
            return _call("DA7", "standard")
        if alt_n < ref_n:
            return _call("DA7", "context")
        return _call("DA7", "auxiliary")
    return None


def _donor_key(pos: int, exon: Exon) -> str | None:
    if exon.strand == "+":
        if pos > exon.end:
            distance = pos - exon.end
            return f"+{distance}" if distance <= 6 else None
        if exon.start < pos <= exon.end:
            distance = exon.end - pos + 1
            return f"E-{distance}" if distance <= 3 else None
        return None
    if pos < exon.start + 1:
        distance = exon.start + 1 - pos
        return f"+{distance}" if distance <= 6 else None
    if exon.start < pos <= exon.end:
        distance = pos - exon.start
        return f"E-{distance}" if distance <= 3 else None
    return None


def _acceptor_key(pos: int, exon: Exon) -> str | None:
    if exon.strand == "+":
        if pos == exon.start + 1:
            return "E+1"
        if pos < exon.start + 1:
            distance = exon.start + 1 - pos
            return f"-{distance}" if distance <= 50 else None
        return None
    if pos == exon.end:
        return "E+1"
    if pos > exon.end:
        distance = pos - exon.end
        return f"-{distance}" if distance <= 50 else None
    return None


def _donor_genomic(offset: int, exon: Exon) -> int:
    if exon.strand == "+":
        return exon.end + offset + 1 if offset < 0 else exon.end + offset
    return exon.start - offset if offset < 0 else exon.start + 1 - offset


def _acceptor_genomic(offset: int, exon: Exon) -> int:
    if exon.strand == "+":
        return exon.start + 1 if offset > 0 else exon.start + 1 + offset
    return exon.end if offset > 0 else exon.end - offset


def _junction_distance(key: str) -> int:
    if key.startswith("E"):
        return int(key[2:]) - 1
    return abs(int(key)) - 1


def _window(fasta: object, chrom: str, coords: dict[int, int]) -> dict[int, str] | None:
    start = min(coords.values())
    end = max(coords.values())
    try:
        seq = fasta.fetch(chrom, start, end)  # type: ignore[attr-defined]
    except (KeyError, ValueError, OSError):
        return None
    if len(seq) != end - start + 1:
        return None
    return {offset: seq[genomic - start] for offset, genomic in coords.items()}


def _hgvs_hint(hgvsc: object) -> tuple[str, str] | None:
    match = _HGVS_OFFSET.search(str(hgvsc or ""))
    if match is None:
        return None
    side = "donor" if match.group("sign") == "+" else "acceptor"
    return side, f"{match.group('sign')}{int(match.group('dist'))}"


def _transcript_id(hgvsc: object) -> str:
    text = str(hgvsc or "").split(":", 1)[0]
    match = _HGVS_TRANSCRIPT.match(text)
    return match.group("id") if match else ""


def _nearest(
    fasta: object,
    exons: list[Exon],
    chrom: str,
    pos: int,
    ref: str,
    hint: tuple[str, str] | None,
) -> tuple[str, str, dict[int, str]] | None:
    best: tuple[tuple[int, int], str, str, dict[int, str]] | None = None
    for exon in exons:
        if exon.chrom != chrom:
            continue
        for side, key_fn, offsets, genomic in (
            ("donor", _donor_key, _DONOR_OFFSETS, _donor_genomic),
            ("acceptor", _acceptor_key, _ACCEPTOR_OFFSETS, _acceptor_genomic),
        ):
            key = key_fn(pos, exon)
            if key is None:
                continue
            coords = {offset: genomic(offset, exon) for offset in offsets}
            bases = _window(fasta, chrom, coords)
            changed = 1 if key == "E+1" else (-int(key[2:]) if key.startswith("E-") else int(key))
            if bases is None or bases.get(changed) != ref:
                continue
            prefer = 0 if hint == (side, key) else 1
            rank = (_junction_distance(key), prefer)
            if best is None or rank < best[0]:
                best = (rank, side, key, bases)
    if best is None:
        return None
    return best[1], best[2], best[3]


def annotate_splice_heuristics(
    variants: pd.DataFrame,
    exons: list[Exon] | None,
    fasta: object | None,
) -> pd.DataFrame:
    """Add the four splice-heuristic columns. Missing inputs leave them blank."""
    frame = _blank(variants)
    if frame.empty or not exons or fasta is None:
        return frame
    by_chrom: dict[str, list[Exon]] = {}
    for exon in exons:
        by_chrom.setdefault(exon.chrom, []).append(exon)

    for index, row in frame.iterrows():
        ref = str(row.get("ref") or "").upper()
        alt = str(row.get("alt") or "").upper()
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
            continue
        pos = row.get("pos")
        if pd.isna(pos):
            continue
        chrom = chrom_bare(row.get("chrom"))
        gene = str(row.get("gene") or "").strip().upper()
        if gene in {"", "NAN", "NONE", "."}:
            gene = ""
        transcript = _transcript_id(row.get("hgvsc"))
        chrom_exons = by_chrom.get(chrom, [])
        transcript_exons = [exon for exon in chrom_exons if transcript and exon.transcript == transcript]
        gene_exons = [exon for exon in chrom_exons if gene and exon.gene == gene]
        hint = _hgvs_hint(row.get("hgvsc"))
        found = _nearest(fasta, transcript_exons or gene_exons, chrom, int(pos), ref, hint)
        if found is None and transcript_exons and gene_exons:
            found = _nearest(fasta, gene_exons, chrom, int(pos), ref, hint)
        if found is None:
            continue
        side, key, bases = found
        called = classify_donor(key, bases, alt) if side == "donor" else classify_acceptor(key, bases, alt)
        if called is None:
            continue
        rule, subgroup = called
        percent, count = _PUBLISHED[(rule, subgroup)]
        frame.at[index, "splice_rule"] = rule
        frame.at[index, "splice_subgroup"] = subgroup
        frame.at[index, "splice_sav_pct"] = percent
        frame.at[index, "splice_heuristic_n"] = count
    return frame


def prepare_exons(path: str | Path, genes: pd.Series | None, aliases: dict[str, str] | None) -> list[Exon]:
    """Load panel exons from a GTF, applying the same gene aliases as the scorer."""
    wanted = None
    if genes is not None:
        wanted = set(apply_gene_aliases(genes, aliases).astype(str))
    return load_gtf_exons(path, wanted, aliases)

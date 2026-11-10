#!/usr/bin/env python3
"""align_conservation_plot.py

Sliding-window conservation plot for pairwise / multiple sequence alignments.

Reads EMBOSS ``srspair`` / ``pair`` output as well as FASTA, Clustal, Stockholm,
PHYLIP, MSF and Nexus alignments, computes per-window identity and gap density,
optionally draws a gene track from GenBank / GFF3 / TSV annotation (or from
automatic ORF prediction), and exports publication-quality SVG / PDF / PNG.

Examples
--------
CLI::

    python align_conservation_plot.py -i aln.srspair -o fig.svg --pdf \
        --coords-seq 2 --annotation NC_004718.3

``--annotation`` takes either a local GenBank/GFF3/TSV file or an NCBI
accession, which is downloaded once and cached under the user cache
directory (``~/Library/Caches/align-conservation-plot`` on macOS).

GUI::

    python align_conservation_plot.py --gui
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import sys
import traceback
import urllib.parse
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------

GAP_CHARS = {"-", "."}

# EMBOSS style:  NAME   START  SEQUENCE  END
SEQ_LINE_RE = re.compile(
    r"^(?P<name>\S.*?)\s+(?P<start>\d+)\s+(?P<seq>[A-Za-z\-\.\*]*)\s+(?P<end>\d+)\s*$"
)


@dataclass
class Alignment:
    names: List[str]
    seqs: List[str]
    source_format: str = "unknown"
    description: str = ""

    @property
    def nseq(self) -> int:
        return len(self.seqs)

    @property
    def ncol(self) -> int:
        return len(self.seqs[0]) if self.seqs else 0

    def ungapped(self, idx: int) -> str:
        return "".join(c for c in self.seqs[idx] if c not in GAP_CHARS)


@dataclass
class WindowStat:
    start: int
    end: int
    ncols: int
    matches: int
    mismatches: int
    gaps: int
    ident_ungapped: float
    ident_with_gaps: float
    gap_frac: float


@dataclass
class IndelStat:
    seq_index: int
    coord: int
    length: int
    kind: str  # "deletion" (gap in this sequence) or "insertion"
    aln_col: int


@dataclass
class Hotspot:
    start: int
    end: int
    ident: float
    gap_frac: float
    ncols: int
    matches: float = 0.0
    mismatches: float = 0.0
    gaps: int = 0
    indel_bp: int = 0
    pairs_per_column: int = 1

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def ident_with_gaps(self) -> float:
        """Identity when gap columns stay in the denominator.

        ``matches`` counts identical sequence pairs, so a column contributes
        ``pairs_per_column`` comparisons.  Without this factor the value is
        inflated by n(n-1)/2 for multi-sequence alignments.
        """
        denom = self.ncols * max(1, self.pairs_per_column)
        return 100.0 * self.matches / denom if denom else 0.0

    @property
    def indel_bp_pct(self) -> float:
        """Share of the region covered by reported (>= threshold) indels."""
        return 100.0 * self.indel_bp / self.length if self.length else 0.0


@dataclass
class Feature:
    name: str
    start: int
    end: int
    strand: str = "+"
    color: Optional[str] = None
    source: str = "annotation"


@dataclass
class AnalysisResult:
    alignment: Alignment
    coords_seq: int
    coord_values: List[int]
    lead_gap: int
    trail_gap: int
    core_columns: int
    total_matches: float
    total_mismatches: float
    total_gaps: int
    indels: List[IndelStat]
    windows: List[WindowStat]
    hotspots: List[Hotspot]
    features: List[Feature]
    feature_source: str
    feature_source_origin: str = ""
    pairs_per_column: int = 1
    warnings: List[str] = field(default_factory=list)

    @property
    def identity_ungapped(self) -> float:
        denom = self.total_matches + self.total_mismatches
        return 100.0 * self.total_matches / denom if denom else 0.0

    @property
    def identity_with_gaps(self) -> float:
        denom = self.core_columns * max(1, self.pairs_per_column)
        return 100.0 * self.total_matches / denom if denom else 0.0

    @property
    def coord_min(self) -> int:
        return min(self.coord_values)

    @property
    def coord_max(self) -> int:
        return max(self.coord_values)


# --------------------------------------------------------------------------
# Alignment parsing
# --------------------------------------------------------------------------


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_emboss_text(text: str) -> Optional[Alignment]:
    """Parse EMBOSS 'srspair' or 'pair' output.

    The format is a sequence of blocks; each block holds coordinate-annotated
    sequence lines separated by a match line (``|``/``:``/``.``).  Names are
    truncated to a fixed width, so they are matched literally across blocks.
    """
    order: List[str] = []
    chunks: Dict[str, List[str]] = {}
    last_end: Dict[str, int] = {}
    saw_header = False

    for raw in text.splitlines():
        if raw.startswith("#"):
            if raw.startswith("# Aligned_sequences") or raw.startswith("# 1:"):
                saw_header = True
            continue
        line = raw.rstrip("\n")
        if not line.strip():
            continue

        m = SEQ_LINE_RE.match(line)
        if m and m.group("seq"):
            name = m.group("name").strip()
            start = int(m.group("start"))
            seq = m.group("seq")
            end = int(m.group("end"))
            if name not in chunks:
                order.append(name)
                chunks[name] = []
                last_end[name] = None
            protein_ok = all(c in "-.*" or c.isalpha() for c in seq)
            if protein_ok:
                prev = last_end[name]
                if prev is not None and start not in (0, prev, prev + 1):
                    return None
                chunks[name].append(seq.upper())
                last_end[name] = end
                continue

        # match / similarity line: no letters or digits at all
        if not re.search(r"[A-Za-z0-9]", line):
            continue

    if not saw_header and len(order) < 2:
        return None
    if len(order) < 2:
        return None

    names = order
    seqs = ["".join(chunks[n]) for n in names]
    if len({len(s) for s in seqs}) != 1 or not seqs[0]:
        return None
    # A valid alignment must contain real residues somewhere.
    if not any(any(c not in GAP_CHARS for c in s) for s in seqs):
        return None
    return Alignment(names=names, seqs=seqs, source_format="emboss")


def parse_with_biopython(path: str) -> Optional[Alignment]:
    try:
        from Bio import AlignIO  # noqa: WPS433 (lazy import)
    except Exception:
        return None

    formats = ["fasta", "clustal", "stockholm", "phylip-relaxed", "phylip", "msf", "nexus"]
    for fmt in formats:
        try:
            with _open_text(path) as handle:
                aln = AlignIO.read(handle, fmt)
            if len(aln) >= 2 and aln.get_alignment_length() > 0:
                return Alignment(
                    names=[r.id for r in aln],
                    seqs=[str(r.seq).upper() for r in aln],
                    source_format=fmt,
                )
        except Exception:
            continue
    return None


def read_alignment(path: str) -> Alignment:
    with _open_text(path) as handle:
        text = handle.read()
    aln = parse_emboss_text(text)
    if aln is not None:
        return aln
    aln = parse_with_biopython(path)
    if aln is None:
        raise ValueError(
            "无法识别的比对格式。支持：EMBOSS srspair/pair、FASTA、Clustal、"
            "Stockholm、PHYLIP、MSF、Nexus。"
        )
    return aln


def read_sequence_names(path: str) -> List[str]:
    """Sequence names exactly as the parser reports them.

    The alias editor must key on these, and srspair truncates names to a fixed
    width, so a header sniff would not match.
    """
    try:
        return read_alignment(path).names
    except Exception:
        return []


# --------------------------------------------------------------------------
# Column statistics
# --------------------------------------------------------------------------


def is_gap(ch: str) -> bool:
    return ch in GAP_CHARS


def column_classes(seqs: Sequence[str]) -> Tuple[List[bool], List[bool], List[int]]:
    """Return (gap_mask_per_seq, all_gap, informative_count)."""
    ncol = len(seqs[0])
    gap_mask = [[is_gap(s[i]) for i in range(ncol)] for s in seqs]
    all_gap = [all(gap_mask[k][i] for k in range(len(seqs))) for i in range(ncol)]
    informative = [sum(1 for k in range(len(seqs)) if not gap_mask[k][i]) for i in range(ncol)]
    return gap_mask, all_gap, informative


def column_pair_score(seqs: Sequence[str], i: int) -> Tuple[float, float]:
    """Return (identical residue pairs, informative residue pairs) in a column.

    For two sequences this is (1, 1) for a match, (0, 1) for a mismatch and
    (0, 0) when a gap is present.  For an MSA it yields the average pairwise
    identity contribution of that column.
    """
    residues = [s[i] for s in seqs if not is_gap(s[i])]
    n = len(residues)
    if n < 2:
        return 0.0, 0.0
    total = n * (n - 1) / 2.0
    counts: Dict[str, int] = {}
    for ch in residues:
        counts[ch] = counts.get(ch, 0) + 1
    same = sum(c * (c - 1) / 2.0 for c in counts.values())
    return same, total


def coord_array(seq: str) -> List[int]:
    """1-based coordinate of every alignment column (gaps inherit previous)."""
    coords: List[int] = []
    pos = 0
    for ch in seq:
        if is_gap(ch):
            coords.append(pos)
        else:
            pos += 1
            coords.append(pos)
    # leading gaps inherit the first real coordinate
    first = next((c for c in coords if c > 0), 0)
    return [c if c > 0 else first for c in coords]


def analyse(
    aln: Alignment,
    coords_seq: int = 1,
    window: int = 300,
    step: int = 50,
    indel_min: int = 8,
    hotspot_k: float = 1.5,
    hotspot_pct: Optional[float] = None,
    keep_terminal_gaps: bool = False,
) -> AnalysisResult:
    coords_seq = coords_seq - 1
    if not 0 <= coords_seq < aln.nseq:
        raise ValueError(f"--coords-seq 超出范围（1..{aln.nseq}）")

    seqs = aln.seqs
    ncol = aln.ncol
    gap_mask, all_gap, informative = column_classes(seqs)

    # Terminal columns where fewer than two sequences carry a residue are not
    # comparable: in a pairwise alignment this is exactly a one-sided overhang
    # (e.g. needle padding the shorter sequence with '-').  Trim them.
    lead = 0
    while lead < ncol and informative[lead] < 2:
        lead += 1
    trail = 0
    while trail < ncol - lead and informative[ncol - 1 - trail] < 2:
        trail += 1

    if keep_terminal_gaps:
        first_col, last_col = 0, ncol - 1
    else:
        first_col, last_col = lead, ncol - 1 - trail
    if last_col < first_col:
        raise ValueError("比对中没有可比对的列（全部为 gap）")

    core_slice = range(first_col, last_col + 1)
    pairs_per_column = aln.nseq * (aln.nseq - 1) // 2

    def column_has_gap(i: int) -> bool:
        return any(gap_mask[k][i] for k in range(aln.nseq))

    match_pairs = 0.0
    pair_total = 0.0
    gaps = 0
    for i in core_slice:
        same, total = column_pair_score(seqs, i)
        if column_has_gap(i):
            gaps += 1
        match_pairs += same
        pair_total += total

    if aln.nseq == 2:
        matches: float = int(round(match_pairs))
        mismatches: float = int(round(pair_total - match_pairs))
    else:
        matches = round(match_pairs, 1)
        mismatches = round(pair_total - match_pairs, 1)

    coords = coord_array(seqs[coords_seq])

    # ---- indels -----------------------------------------------------------
    indels: List[IndelStat] = []
    for k in range(aln.nseq):
        i = first_col
        while i <= last_col:
            if gap_mask[k][i]:
                j = i
                while j + 1 <= last_col and gap_mask[k][j + 1]:
                    j += 1
                length = j - i + 1
                if length >= indel_min:
                    others = any(not gap_mask[o][i] for o in range(aln.nseq) if o != k)
                    if others:
                        anchor = coords[i - 1] if i > 0 else coords[i]
                        indels.append(
                            IndelStat(
                                seq_index=k,
                                coord=anchor,
                                length=length,
                                kind="deletion",
                                aln_col=i,
                            )
                        )
                i = j + 1
            else:
                i += 1

    # ---- windows over reference coordinates -------------------------------
    coord_min = min(coords[i] for i in core_slice)
    coord_max = max(coords[i] for i in core_slice)
    windows: List[WindowStat] = []
    by_coord: Dict[int, List[int]] = {}
    for i in core_slice:
        by_coord.setdefault(coords[i], []).append(i)

    start = coord_min
    last_start = max(coord_min, coord_max - window + 1)
    starts = list(range(coord_min, last_start + 1, step)) or [coord_min]
    for start in starts:
        end = start + window - 1
        cols: List[int] = []
        for c in range(start, end + 1):
            cols.extend(by_coord.get(c, ()))
        if not cols:
            continue
        wm = 0.0
        wtot = 0.0
        wg = 0
        for i in cols:
            same, total = column_pair_score(seqs, i)
            if column_has_gap(i):
                wg += 1
            wm += same
            wtot += total
        windows.append(
            WindowStat(
                start=start,
                end=end,
                ncols=len(cols),
                matches=int(round(wm)) if aln.nseq == 2 else round(wm, 1),
                mismatches=int(round(wtot - wm)) if aln.nseq == 2 else round(wtot - wm, 1),
                gaps=wg,
                ident_ungapped=100.0 * wm / wtot if wtot else float("nan"),
                ident_with_gaps=100.0 * wm / (len(cols) * max(1, pairs_per_column)),
                gap_frac=100.0 * wg / len(cols),
            )
        )

    # ---- hotspots ---------------------------------------------------------
    usable = [
        w
        for w in windows
        if not math.isnan(w.ident_ungapped) and w.ncols >= max(20, window * 0.5)
    ]
    hotspots: List[Hotspot] = []
    if len(usable) >= 3:
        values = sorted(w.ident_ungapped for w in usable)
        if hotspot_pct is not None:
            idx = max(0, min(len(values) - 1, int(round(hotspot_pct / 100.0 * (len(values) - 1)))))
            threshold = values[idx]
        else:
            med = statistics_median(values)
            mad = statistics_median([abs(v - med) for v in values]) or 0.5
            threshold = med - hotspot_k * 1.4826 * mad
        low = [w for w in usable if w.ident_ungapped <= threshold]
        merged: List[List[WindowStat]] = []
        for w in low:
            if merged and w.start <= merged[-1][-1].end + 1:
                merged[-1].append(w)
            else:
                merged.append([w])
        def region_stats(start: int, end: int) -> Tuple[float, float, int, int]:
            """Exact column statistics for a coordinate range (windows overlap)."""
            cols: List[int] = []
            for c in range(start, end + 1):
                cols.extend(by_coord.get(c, ()))
            same_total = pair_total = 0.0
            gap_cols = 0
            for i in cols:
                same, total = column_pair_score(seqs, i)
                if column_has_gap(i):
                    gap_cols += 1
                same_total += same
                pair_total += total
            return same_total, pair_total - same_total, gap_cols, len(cols)

        for group in merged:
            if len(group) < 2 and len(low) > 1:
                # keep isolated single windows only when they are clearly worst
                if group[0].ident_ungapped > threshold - 1.0:
                    continue
            hs_start, hs_end = group[0].start, group[-1].end
            m, mm, gp, n = region_stats(hs_start, hs_end)
            indel_bp = sum(
                ind.length for ind in indels if hs_start <= ind.coord <= hs_end
            )
            hotspots.append(
                Hotspot(
                    start=hs_start,
                    end=hs_end,
                    ident=100.0 * m / (m + mm) if m + mm else 0.0,
                    gap_frac=100.0 * gp / n if n else 0.0,
                    ncols=n,
                    matches=m,
                    mismatches=mm,
                    gaps=gp,
                    indel_bp=indel_bp,
                    pairs_per_column=pairs_per_column,
                )
            )

    return AnalysisResult(
        alignment=aln,
        coords_seq=coords_seq,
        coord_values=coords[first_col : last_col + 1],
        lead_gap=lead,
        trail_gap=trail,
        core_columns=last_col - first_col + 1,
        total_matches=matches,
        total_mismatches=mismatches,
        total_gaps=gaps,
        indels=indels,
        windows=windows,
        hotspots=hotspots,
        features=[],
        feature_source="none",
        pairs_per_column=pairs_per_column,
    )


def statistics_median(values: Sequence[float]) -> float:
    s = sorted(values)
    n = len(s)
    if n == 0:
        return 0.0
    if n % 2:
        return s[n // 2]
    return 0.5 * (s[n // 2 - 1] + s[n // 2])


# --------------------------------------------------------------------------
# Annotation
# --------------------------------------------------------------------------

PALETTE = [
    "#cfe3f7", "#b9d4ef", "#f6b3a6", "#f8cfa0", "#c8e6c9",
    "#a5d6a7", "#dcedc8", "#b2dfdb", "#80cbc4", "#ffcc80",
    "#ce93d8", "#e0e0e0", "#b0bec5",
]


def load_gff3(path: str) -> List[Feature]:
    feats: List[Feature] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype, start, end, strand, attrs = parts[2], parts[3], parts[4], parts[6], parts[8]
            if ftype.lower() not in {"gene", "cds", "trna", "rrna", "misc_feature"}:
                continue
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            name = None
            for key in ("Name", "gene", "locus_tag", "ID", "Parent"):
                m = re.search(key + r"=([^;]+)", attrs)
                if m:
                    name = m.group(1)
                    break
            feats.append(Feature(name or f"{ftype}:{s}-{e}", s, e, strand if strand in "+-" else "+"))
    return feats


def load_tsv(path: str) -> List[Feature]:
    feats: List[Feature] = []
    with _open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 3:
                parts = [p.strip() for p in re.split(r"[,\s]+", line) if p.strip()]
            if len(parts) < 3:
                continue
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            strand = parts[3] if len(parts) > 3 and parts[3] in "+-" else "+"
            color = parts[4] if len(parts) > 4 and parts[4].startswith("#") else None
            feats.append(Feature(parts[0], s, e, strand, color))
    return feats


def _dedupe_features(feats: List[Feature]) -> List[Feature]:
    seen = set()
    out: List[Feature] = []
    for f in sorted(feats, key=lambda x: (x.start, x.end)):
        key = (f.name, f.start, f.end, f.strand)
        if key in seen:
            continue
        seen.add(key)
        out.append(f)
    return out


def _name_from_qualifiers(qualifiers: Dict[str, List[str]], fallback: str) -> str:
    for key in ("gene", "locus_tag", "product", "note"):
        values = qualifiers.get(key)
        if values:
            value = str(values[0]).strip()
            if value:
                return value
    return fallback


def load_genbank(path: str, record_id: Optional[str] = None) -> List[Feature]:
    handle_ctx = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")
    with handle_ctx as handle:
        return parse_genbank_text(handle.read(), record_id, source_name=path)


def parse_genbank_text(
    text: str, record_id: Optional[str] = None, source_name: str = "<genbank>"
) -> List[Feature]:
    import io

    from Bio import SeqIO

    records = list(SeqIO.parse(io.StringIO(text), "genbank"))
    if not records:
        raise ValueError(f"{source_name} 不是有效的 GenBank 记录")
    record = records[0]
    if record_id:
        for r in records:
            if record_id in (r.id, r.name, r.description):
                record = r
                break

    buckets: Dict[str, List[Feature]] = {"gene": [], "CDS": [], "other": []}
    for f in record.features:
        if f.type in {"gene"}:
            bucket = "gene"
        elif f.type in {"CDS", "tRNA", "rRNA"}:
            bucket = "CDS"
        elif f.type in {"mat_peptide", "misc_feature"}:
            bucket = "other"
        else:
            continue
        try:
            s = int(f.location.start) + 1
            e = int(f.location.end)
        except Exception:
            continue
        strand = "-" if f.location.strand == -1 else "+"
        name = _name_from_qualifiers(f.qualifiers, f.type)
        buckets[bucket].append(Feature(name, s, e, strand))

    for key in ("gene", "CDS", "other"):
        if buckets[key]:
            return _dedupe_features(buckets[key])
    return []


ACCESSION_RE = re.compile(r"^[A-Za-z]{1,6}_?\d+(?:\.\d+)?$")
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class AnnotationFetchError(RuntimeError):
    """Raised when an online annotation cannot be retrieved."""


def looks_like_accession(value: str) -> bool:
    """True for NCBI-style accessions such as NC_004718.3 or MN908947."""
    value = (value or "").strip()
    if not value or os.path.exists(value):
        return False
    if "://" in value or os.sep in value:
        return False
    return bool(ACCESSION_RE.match(value))


def annotation_cache_dir() -> str:
    override = os.environ.get("ALIGN_PLOT_CACHE")
    if override:
        return os.path.expanduser(override)
    if sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Caches")
    elif os.name == "nt":
        base = os.environ.get("LOCALAPPDATA", os.path.expanduser("~"))
    else:
        base = os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
    return os.path.join(base, "align-conservation-plot")


def cache_path_for(accession: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", accession.upper())
    return os.path.join(annotation_cache_dir(), f"{safe}.gb")


def _clean_ncbi_detail(raw: str) -> str:
    """Turn NCBI's percent/plus-encoded error body into readable text."""
    text = urllib.parse.unquote(raw)
    text = re.sub(r"(?<=\w)\+(?=\w)", "", text)
    text = text.replace("+", " ")
    text = " ".join(text.split())
    text = re.sub(r"(?<=\d)\s*\.\s*(?=\d)", ".", text)
    text = re.sub(r"(?<=\w)\s+:\s*(?=\w)", ": ", text)
    return text[:200]


def _http_get(url: str, timeout: float) -> str:
    import urllib.error
    import urllib.request

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "align-conservation-plot/1.1 (+https://github.com/)",
            "Accept": "text/plain",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        charset = response.headers.get_content_charset() or "utf-8"
        return response.read().decode(charset, errors="replace")


def fetch_genbank_text(
    accession: str,
    timeout: float = 60.0,
    use_cache: bool = True,
    refresh: bool = False,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
    retries: int = 3,
) -> Tuple[str, str]:
    """Download a GenBank record from NCBI.

    Returns ``(genbank_text, origin)`` where origin is ``"cache"`` or
    ``"download"``.  Network access honours ``http_proxy`` / ``https_proxy``
    from the environment, so a local proxy is honoured automatically.
    """
    import time
    import urllib.error

    accession = accession.strip().upper()
    path = cache_path_for(accession)
    if use_cache and not refresh and os.path.exists(path):
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
        if "LOCUS" in text:
            return text, "cache"

    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gbwithparts",
        "retmode": "text",
    }
    email = email or os.environ.get("NCBI_EMAIL")
    api_key = api_key or os.environ.get("NCBI_API_KEY")
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    query = "&".join(f"{k}={urllib.parse.quote(str(v))}" for k, v in params.items())
    url = f"{NCBI_EFETCH}?{query}"

    last_error: Optional[Exception] = None
    for attempt in range(retries):
        try:
            text = _http_get(url, timeout)
        except urllib.error.HTTPError as exc:
            detail = ""
            try:
                detail = exc.read().decode("utf-8", errors="replace")
            except Exception:
                pass
            detail = _clean_ncbi_detail(detail)
            if exc.code in (400, 404) or "Failed to retrieve" in detail:
                raise AnnotationFetchError(
                    f"NCBI 找不到 accession {accession}（HTTP {exc.code}）。"
                    "请确认编号与版本号正确，例如 NC_004718.3。"
                    + (f" 服务器返回：{detail}" if detail else "")
                ) from exc
            raise AnnotationFetchError(
                f"NCBI 返回 HTTP {exc.code}（accession {accession}）。{detail}".strip()
            ) from exc
        except Exception as exc:  # noqa: BLE001 - URLError, timeout, DNS, proxy
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(1.5 * (attempt + 1))
                continue
            break

        if "LOCUS" in text and "FEATURES" in text:
            if use_cache:
                try:
                    os.makedirs(os.path.dirname(path), exist_ok=True)
                    with open(path, "w", encoding="utf-8") as handle:
                        handle.write(text)
                except OSError:
                    pass
            return text, "download"

        compact = " ".join(text.split())[:300]
        if "Error" in text or "error" in text or not compact:
            raise AnnotationFetchError(
                f"NCBI 未返回 GenBank 记录（accession {accession}）：{compact or '空响应'}"
            )
        raise AnnotationFetchError(
            f"NCBI 响应中没有 FEATURES 段（accession {accession}）：{compact}"
        )

    raise AnnotationFetchError(
        f"无法连接 NCBI 获取 {accession}：{last_error}。"
        "请检查网络；若本机通过代理上网，可设置 "
        "http_proxy/https_proxy（例如 http://127.0.0.1:8080）。"
    )


def load_genbank_accession(
    accession: str,
    record_id: Optional[str] = None,
    timeout: float = 60.0,
    use_cache: bool = True,
    refresh: bool = False,
) -> Tuple[List[Feature], str]:
    text, origin = fetch_genbank_text(
        accession, timeout=timeout, use_cache=use_cache, refresh=refresh
    )
    feats = parse_genbank_text(text, record_id, source_name=accession)
    return feats, origin


STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _six_frame_orfs(seq: str, min_aa: int) -> List[Tuple[int, int, str]]:
    """Return (start, end, strand) 1-based inclusive ORFs starting with ATG."""
    seq = seq.upper()
    rc = seq.translate(_COMPLEMENT)[::-1]
    out: List[Tuple[int, int, str]] = []
    for strand, s in (("+", seq), ("-", rc)):
        for frame in range(3):
            i = frame
            while i + 3 <= len(s):
                if s[i : i + 3] != "ATG":
                    i += 3
                    continue
                j = i
                while j + 3 <= len(s):
                    if STANDARD_CODE.get(s[j : j + 3]) == "*":
                        break
                    j += 3
                n_aa = (j - i) // 3
                if n_aa >= min_aa:
                    if strand == "+":
                        out.append((i + 1, j + 2, "+"))
                    else:
                        out.append((len(s) - j + 1, len(s) - i, "-"))
                i = j + 3
    return out


def predict_orfs(seq: str, min_aa: int = 90) -> List[Feature]:
    orfs = _six_frame_orfs(seq, min_aa)
    orfs.sort(key=lambda t: (-(t[1] - t[0]), t[0]))
    kept: List[Tuple[int, int, str]] = []
    for s, e, strand in orfs:
        if all(e < ks or s > ke for ks, ke, _ in kept):
            kept.append((s, e, strand))
    kept.sort()
    return [Feature(f"ORF{i + 1}", s, e, strand, source="predicted") for i, (s, e, strand) in enumerate(kept)]


def attach_features(
    result: AnalysisResult,
    path: Optional[str],
    fmt: Optional[str] = None,
    record_id: Optional[str] = None,
    min_orf_aa: int = 90,
    annotation_options: Optional[dict] = None,
) -> None:
    aln = result.alignment
    if path:
        if looks_like_accession(path):
            opts = annotation_options or {}
            if fmt:
                result.warnings.append("--annotation-format is ignored for accession input")
            feats, origin = load_genbank_accession(
                path,
                record_id,
                timeout=opts.get("timeout", 60.0),
                use_cache=opts.get("use_cache", True),
                refresh=opts.get("refresh", False),
            )
            result.features = feats
            result.feature_source = f"GenBank {path.strip().upper()}"
            result.feature_source_origin = origin
            if not feats:
                result.warnings.append(f"no usable features in {path}")
            return

        if not os.path.exists(path):
            raise FileNotFoundError(2, "annotation file not found", path)

        ext = os.path.splitext(path)[1].lower()
        fmt = (fmt or ("genbank" if ext in (".gb", ".gbk", ".genbank", ".gbff") else
                       "gff3" if ext in (".gff", ".gff3") else "tsv")).lower()
        if fmt in {"gb", "gbk", "genbank"}:
            feats = load_genbank(path, record_id)
            result.feature_source = "GenBank"
        elif fmt in {"gff", "gff3"}:
            feats = load_gff3(path)
            result.feature_source = "GFF3"
        else:
            feats = load_tsv(path)
            result.feature_source = "TSV"
        result.features = feats
        if not feats:
            result.warnings.append(
                f"no usable features in {os.path.basename(path)}"
            )
        return

    seq = aln.ungapped(result.coords_seq)
    if len(seq) >= 300 and set(seq.upper()) <= set("ACGTN"):
        result.features = predict_orfs(seq, min_orf_aa)
        result.feature_source = "predicted ORFs"
    else:
        result.feature_source = "none"


def feature_source_label(result: AnalysisResult, lang: str = "en") -> str:
    """Render the annotation source, including how an accession was obtained."""
    source = result.feature_source
    origin = getattr(result, "feature_source_origin", "")
    if not origin:
        return source
    if lang == "zh":
        suffix = "已下载" if origin == "download" else "缓存"
    else:
        suffix = "downloaded" if origin == "download" else "cached"
    return f"{source} ({suffix})"


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

LANG_PRESETS = {
    "en": {
        "identity": "Identity (%)",
        "gap": "Gap (%)",
        "pos": "Position in {name} (bp)",
        "genes": "genes",
        "median": "median {value:.1f}%",
        "legend_identity": "identity (gap-free)",
        "legend_hotspot": "low-conservation hotspot",
        "subtitle_columns": "{value:,} alignment columns",
        "subtitle_identity": "{value:.1f}% identity",
        "subtitle_multi": "{nseq} sequences, mean pairwise identity",
        "footer": (
            "Window {window} bp / step {step} bp. Identity excludes gap columns; "
            "median {median:.1f}%. Terminal all-gap columns {lead:,} / {trail:,} "
            "were excluded. Genes: {source}. Compared region: {start:,}-{end:,} "
            "({span:,} bp)."
        ),
    },
    "zh": {
        "identity": "一致度 (%)",
        "gap": "gap 占比 (%)",
        "pos": "{name} 上的位置 (bp)",
        "genes": "基因",
        "median": "中位数 {value:.1f}%",
        "legend_identity": "一致度（不计 gap）",
        "legend_hotspot": "低保守热点",
        "subtitle_columns": "{value:,} 列",
        "subtitle_identity": "{value:.1f}% 一致度",
        "subtitle_multi": "{nseq} 条序列，平均两两一致度",
        "footer": (
            "窗口 {window} bp / 步长 {step} bp；一致度不计 gap 列，中位数 {median:.1f}%。"
            "已剔除两端全 gap 列 {lead:,} / {trail:,} 列。基因来源：{source}。"
            "比较区间：{start:,}–{end:,}（{span:,} bp）。"
        ),
    },
}


# Font choices, ordered by three rules: open licence first, sans-serif before
# serif, CJK-capable before Latin-only.  See FONTS.md for the licence evidence.
FONT_LICENSES: Dict[str, Tuple[str, bool]] = {
    # SIL Open Font License 1.1 — embeddable and redistributable
    "Source Han Sans SC": ("SIL Open Font License 1.1", True),
    "Source Han Sans": ("SIL Open Font License 1.1", True),
    "Noto Sans CJK SC": ("SIL Open Font License 1.1", True),
    "Noto Sans CJK": ("SIL Open Font License 1.1", True),
    "Noto Sans": ("SIL Open Font License 1.1", True),
    "Source Sans Pro": ("SIL Open Font License 1.1", True),
    "Source Sans 3": ("SIL Open Font License 1.1", True),
    "Source Han Serif SC": ("SIL Open Font License 1.1", True),
    # Bitstream Vera licence — permissive, embeddable
    "DejaVu Sans": ("Bitstream Vera licence", True),
    # Proprietary system fonts: usable on the machine that has them, but they
    # cannot be redistributed with the project.
    "PingFang SC": ("proprietary licence (Apple/DynaComware)", False),
    "PingFang HK": ("proprietary licence (Apple/DynaComware)", False),
    "PingFang TC": ("proprietary licence (Apple/DynaComware)", False),
    "Hiragino Sans GB": ("proprietary licence (Apple/SCREEN)", False),
    "Heiti SC": ("proprietary licence (Apple/SinoType)", False),
    "Microsoft YaHei": ("proprietary licence (Microsoft)", False),
    "SimHei": ("proprietary licence (Microsoft)", False),
    "Arial Unicode MS": ("proprietary licence (Monotype)", False),
    "Songti SC": ("proprietary licence (Apple/SinoType)", False),
}

# Render order: open + sans + CJK, then open + sans, then proprietary sans,
# then proprietary serif as the final resort.
FIGURE_FONT_ORDER: Tuple[str, ...] = (
    "Source Han Sans SC",
    "Source Han Sans",
    "Noto Sans CJK SC",
    "Noto Sans CJK",
    "Noto Sans",
    "Source Sans Pro",
    "Source Sans 3",
    "DejaVu Sans",
    "PingFang SC",
    "PingFang HK",
    "PingFang TC",
    "Hiragino Sans GB",
    "Heiti SC",
    "Microsoft YaHei",
    "SimHei",
    "Arial Unicode MS",
    "Songti SC",
    "Source Han Serif SC",
)

UI_FONT_CANDIDATES = FIGURE_FONT_ORDER


def pick_ui_font(available: Iterable[str]) -> str:
    """First installed UI font from the cross-platform candidate list."""
    installed = set(available)
    for family in UI_FONT_CANDIDATES:
        if family in installed:
            return family
    return "TkDefaultFont"


def pick_figure_font(available: Iterable[str]) -> Optional[str]:
    """First installed figure font, in licence/sans/CJK preference order."""
    installed = set(available)
    for family in FIGURE_FONT_ORDER:
        if family in installed:
            return family
    return None


def font_license(name: str) -> Tuple[str, bool]:
    return FONT_LICENSES.get(name, ("unknown licence", False))


def font_notice(name: Optional[str], lang: str = "en") -> str:
    """One sentence naming the figure font and its licence, for CLI and GUI."""
    if not name:
        if lang == "zh":
            return "未找到候选字体，将回退到 matplotlib 自带字体；许可信息未知，请自行核实。"
        return (
            "No preferred font was found, so matplotlib's built-in font will be used; "
            "its licence is unknown here, so please verify it yourself."
        )
    licence, open_licence = font_license(name)
    if lang == "zh":
        kind = "开源" if open_licence else "专有"
        return (
            f"图表将使用「{name}」（{licence}，{kind}许可）渲染。"
            "「数据」与「低保守区表 Markdown」导出不含渲染文字。"
            "许可信息可能已过时，正式发布前请自行核实。"
        )
    kind = "open licence" if open_licence else "proprietary licence"
    return (
        f'Figures will be rendered with "{name}" ({licence}, {kind}).'
        ' The "Data" and "Hotspot table Markdown" exports contain no rendered text.'
        " Licence information may be outdated, so please verify it before publishing."
    )


def available_font_families() -> List[str]:
    """Family names matplotlib can currently see."""
    from matplotlib import font_manager

    return sorted({f.name for f in font_manager.fontManager.ttflist})


def configure_matplotlib(lang: str) -> None:
    import logging

    import matplotlib

    # PingFang SC has no separate bold face; matplotlib falls back silently but
    # logs a noisy warning on every run.
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

    # A list in font.family is resolved per glyph; font.sans-serif would pick a
    # single font for the whole string and drop characters it lacks.
    matplotlib.rcParams["font.family"] = list(FIGURE_FONT_ORDER)
    matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["svg.fonttype"] = "none"  # keep text editable
    matplotlib.rcParams["pdf.fonttype"] = 42


def _nice_ylim(values: Iterable[float], pad: float = 3.0) -> Tuple[float, float]:
    vals = [v for v in values if not math.isnan(v)]
    if not vals:
        return 0.0, 100.0
    lo = max(0.0, min(vals) - pad)
    hi = min(100.0, max(vals) + pad)
    if hi - lo < 10:
        lo = max(0.0, lo - 5)
        hi = min(100.0, hi + 5)
    return lo, hi


def build_figure(result: AnalysisResult, opts) -> "object":
    from matplotlib.figure import Figure
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter

    t = LANG_PRESETS.get(opts.lang, LANG_PRESETS["en"])
    aln = result.alignment
    names = display_names(result, opts)
    ref_name = names[result.coords_seq]
    xmin, xmax = result.coord_min, result.coord_max

    show_title = not getattr(opts, "no_title", False)
    show_subtitle = not getattr(opts, "no_subtitle", False)
    show_head = show_title or show_subtitle
    show_gap = not getattr(opts, "no_gap_track", False)
    show_gene = not getattr(opts, "no_gene_track", False) and bool(result.features)
    show_footer = not getattr(opts, "no_footer", False)

    # The layout is assembled from whichever tracks are switched on.
    layout: List[Tuple[str, float]] = []
    if show_head:
        layout.append(("head", 0.62))
    layout.append(("identity", 3.4))
    if show_gap:
        layout.append(("gap", 0.85))
    if show_gene:
        layout.append(("gene", 0.9))

    if opts.figsize:
        figsize = tuple(opts.figsize)
    else:
        extra = (0.55 if show_head else 0.0) + (0.85 if show_gap else 0.0) + (0.9 if show_gene else 0.0)
        figsize = (11.0, 3.5 + extra)
    fig = Figure(figsize=figsize, dpi=opts.dpi, layout="constrained")

    gs = fig.add_gridspec(
        nrows=len(layout),
        ncols=1,
        height_ratios=[h for _, h in layout],
        hspace=0.10,
    )
    slots = {name: gs[i] for i, (name, _) in enumerate(layout)}
    ax_head = None
    if show_head:
        ax_head = fig.add_subplot(slots["head"])
        ax_head.axis("off")
    ax = fig.add_subplot(slots["identity"], sharex=ax_head)
    ax_gap = fig.add_subplot(slots["gap"], sharex=ax) if show_gap else None

    # --- hotspot shading ---------------------------------------------------
    for hs in result.hotspots:
        ax.axvspan(hs.start, hs.end, color="#fdd", zorder=0)

    # --- identity curve ----------------------------------------------------
    wins = [w for w in result.windows if not math.isnan(w.ident_ungapped)]
    xs = [(w.start + w.end) / 2 for w in wins]
    ys = [w.ident_ungapped for w in wins]
    ax.plot(xs, ys, color="#1a6fb5", lw=1.6, zorder=3)
    if ys:
        ax.axhline(result.identity_ungapped, color="#888", ls=(0, (4, 3)), lw=0.9, zorder=2)
    ax.set_ylabel(t["identity"], fontsize=10)
    ylo, yhi = _nice_ylim(ys)
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(xmin, xmax)
    ax.grid(axis="y", color="#e8e8e8", lw=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # hotspot numbers: every hotspot is numbered so the table indices match
    for i, hs in enumerate(result.hotspots, start=1):
        offset = 0.06 if (i % 2 or len(result.hotspots) <= 8) else 0.15
        ax.text(
            (hs.start + hs.end) / 2,
            yhi - (yhi - ylo) * offset,
            str(i),
            ha="center",
            va="top",
            fontsize=7.5,
            color="#c0392b",
        )

    handles = [
        Line2D([0], [0], color="#1a6fb5", lw=1.6, label=t["legend_identity"]),
    ]
    labels = [t["legend_identity"]]
    if result.hotspots:
        handles.append(Rectangle((0, 0), 1, 1, color="#fdd", ec="#e8a0a0"))
        labels.append(t["legend_hotspot"])
    ax.legend(
        handles,
        labels,
        loc="upper left",
        fontsize=7.5,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=0.85,
        ncol=2,
        handlelength=1.6,
        borderaxespad=0.3,
    )
    ax.text(
        xmax,
        result.identity_ungapped + (yhi - ylo) * 0.02,
        t["median"].format(value=result.identity_ungapped),
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#666",
    )

    # --- gap track ---------------------------------------------------------
    if ax_gap is not None:
        gap_x = [(w.start + w.end) / 2 for w in wins]
        gap_y = [w.gap_frac for w in wins]
        ax_gap.fill_between(gap_x, gap_y, step="mid", color="#d95f02", alpha=0.35, lw=0)
        ax_gap.plot(gap_x, gap_y, drawstyle="steps-mid", color="#d95f02", lw=0.9)
        ax_gap.set_ylabel(t["gap"], fontsize=9)
        ax_gap.set_ylim(0, max(5.0, (max(gap_y) if gap_y else 0) * 1.25))
        ax_gap.grid(axis="y", color="#eee", lw=0.6)
        ax_gap.set_axisbelow(True)
        for side in ("top", "right"):
            ax_gap.spines[side].set_visible(False)

        # indel ticks belong to the gap track and disappear with it
        if result.indels:
            ymax = ax_gap.get_ylim()[1]
            for ind in result.indels:
                ax_gap.plot([ind.coord, ind.coord], [0, ymax * 0.18], color="#6a3d9a", lw=1.0, alpha=0.8)

    # --- gene track --------------------------------------------------------
    ax_gene = None
    if show_gene:
        ax_gene = fig.add_subplot(slots["gene"], sharex=ax)
        ax_gene.set_ylim(-1.05, 1.05)
        ax_gene.set_yticks([])
        for side in ("top", "right", "left"):
            ax_gene.spines[side].set_visible(False)
        ax_gene.set_ylabel(t["genes"], fontsize=9, rotation=0, ha="right", va="center", labelpad=14)

        span = max(1, xmax - xmin)
        fig_width_pt = figsize[0] * 72.0
        draws = []
        for i, feat in enumerate(result.features):
            s = max(xmin, feat.start)
            e = min(xmin + span, feat.end)
            if e <= s:
                continue
            draws.append((s, e, feat, feat.color or PALETTE[i % len(PALETTE)]))
        draws.sort(key=lambda x: (x[0], x[1]))

        # stack overlapping features into rows (interval graph colouring)
        rows_end: List[float] = []
        placed = []
        for s, e, feat, color in draws:
            for r, last_end in enumerate(rows_end):
                if s > last_end:
                    rows_end[r] = e
                    placed.append((s, e, feat, color, r))
                    break
            else:
                rows_end.append(e)
                placed.append((s, e, feat, color, len(rows_end) - 1))

        n_rows = max(1, len(rows_end))
        row_h = min(0.24, 0.95 / n_rows * 0.55)
        gap = (1.0 - n_rows * row_h) / max(1, n_rows)
        legend_items: List[str] = []
        counter = 0
        for s, e, feat, color, r in placed:
            y0 = 1.0 - (r + 1) * row_h - r * gap
            ax_gene.add_patch(Rectangle((s, y0), e - s, row_h, facecolor=color, edgecolor="#666", lw=0.4))
            tip = min(max((e - s) * 0.45, 2), span * 0.004)
            if feat.strand == "-":
                ax_gene.add_patch(
                    Polygon([(s, y0 + row_h / 2), (s + tip, y0), (s + tip, y0 + row_h)],
                            closed=True, facecolor=color, edgecolor="#666", lw=0.4)
                )
            else:
                ax_gene.add_patch(
                    Polygon([(e, y0 + row_h / 2), (e - tip, y0), (e - tip, y0 + row_h)],
                            closed=True, facecolor=color, edgecolor="#666", lw=0.4)
                )
            width_pt = (e - s) / span * fig_width_pt
            est_pt = len(feat.name) * 4.1 + 5
            if width_pt > est_pt * 1.15:
                ax_gene.text((s + e) / 2, y0 + row_h / 2, feat.name,
                             ha="center", va="center", fontsize=7.2)
            else:
                counter += 1
                ax_gene.text((s + e) / 2, y0 + row_h / 2, str(counter),
                             ha="center", va="center", fontsize=5.8, color="#7b241c")
                legend_items.append(f"{counter} {feat.name}")

        if legend_items:
            lines: List[str] = []
            current = ""
            for item in legend_items:
                candidate = (current + "   " + item).strip()
                if len(candidate) > 118 and current:
                    lines.append(current)
                    current = item
                else:
                    current = candidate
            if current:
                lines.append(current)
            for i, line in enumerate(lines[:3]):
                ax_gene.text(xmin, 0.02 - i * 0.30, line,
                             fontsize=6.2, color="#444", ha="left", va="top")

    # --- axes labels -------------------------------------------------------
    bottom_ax = ax_gene if ax_gene is not None else (ax_gap if ax_gap is not None else ax)
    bottom_ax.set_xlabel(t["pos"].format(name=ref_name), fontsize=9.5)

    def fmt_coord(v, _pos):
        return f"{int(v):,}"

    bottom_ax.xaxis.set_major_formatter(FuncFormatter(fmt_coord))

    # --- titles and footer -------------------------------------------------
    if show_head:
        if show_title:
            if opts.title:
                title = opts.title
            elif aln.nseq == 2:
                title = f"{names[0]} vs {names[1]}"
            else:
                title = f"{names[0]} (+{aln.nseq - 1} more)"
            ax_head.text(0.0, 0.92, title, fontsize=12, fontweight="bold", ha="left", va="top")
        if show_subtitle:
            subtitle_bits = []
            if aln.nseq > 2:
                subtitle_bits.append(t["subtitle_multi"].format(nseq=aln.nseq))
            subtitle_bits.append(t["subtitle_columns"].format(value=result.core_columns))
            subtitle_bits.append(t["subtitle_identity"].format(value=result.identity_ungapped))
            y = 0.18 if show_title else 0.55
            ax_head.text(0.0, y, "  |  ".join(subtitle_bits), fontsize=7.5, color="#555",
                         ha="left", va="bottom")

    if show_footer:
        footer = t["footer"].format(
            window=opts.window,
            step=opts.step,
            median=result.identity_ungapped,
            lead=result.lead_gap,
            trail=result.trail_gap,
            source=feature_source_label(result, opts.lang) if show_gene else "—",
            start=xmin,
            end=xmax,
            span=xmax - xmin + 1,
        )
        fig.text(0.01, -0.02, footer, fontsize=6.5, color="#666", ha="left", va="top")
    return fig


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

HOTSPOT_TABLE_COLUMNS = (
    "index",
    "coords",
    "length",
    "ident_ungapped",
    "ident_with_gaps",
    "indel",
    "region",
)

HOTSPOT_TABLE_HEADERS = {
    "zh": (
        "序号",
        "坐标 (bp)",
        "长度 (bp)",
        "一致度 % (不计gap)",
        "一致度 %（计gap）",
        "插入或删除 % (gap 列 / 插入或删除 bp)",
        "所在区域",
    ),
    "en": (
        "#",
        "Coordinates (bp)",
        "Length (bp)",
        "Identity % (ungapped)",
        "Identity % (with gaps)",
        "Indel % (gap cols / indel bp)",
        "Region",
    ),
}

HOTSPOT_TABLE_NOTE = {
    "zh": (
        "插入或删除 % 同时给出两个数：区间内 gap 列占比 / "
        "该区间内 ≥ 阈值的插入或删除事件覆盖的比例。"
    ),
    "en": (
        "Indel % shows gap-column share / share covered by indels at or above "
        "the reporting threshold."
    ),
}


# Columns a user may hide; "#" and coordinates always stay.
HIDEABLE_TABLE_COLUMNS = ("length", "ident_ungapped", "ident_with_gaps", "indel", "region")
ALWAYS_VISIBLE_TABLE_COLUMNS = ("index", "coords")


def visible_table_columns(opts) -> Tuple[str, ...]:
    """Table columns left after the display switches, in canonical order."""
    hidden = set(getattr(opts, "hide_table_column", None) or ())
    if getattr(opts, "no_table_note", False):
        # The note explains the indel figure, so they are shown or hidden together.
        hidden.add("indel")
    if getattr(opts, "no_gene_track", False):
        hidden.add("region")
    return tuple(
        col
        for col in HOTSPOT_TABLE_COLUMNS
        if col in ALWAYS_VISIBLE_TABLE_COLUMNS or col not in hidden
    )


def display_options(
    *,
    fig_title: bool = True,
    fig_subtitle: bool = True,
    gap_track: bool = True,
    gene_track: bool = True,
    fig_footer: bool = True,
    table_title: bool = True,
    table_subtitle: bool = True,
    table_note: bool = True,
    col_length: bool = True,
    col_ident_ungapped: bool = True,
    col_ident_gaps: bool = True,
    col_indel: bool = True,
    col_region: bool = True,
) -> Dict[str, object]:
    """Map positive "show it" switches onto the negative CLI flags.

    Keeping this pure lets the GUI build parser-compatible options without
    duplicating the flag names.
    """
    hidden = []
    if not col_length:
        hidden.append("length")
    if not col_ident_ungapped:
        hidden.append("ident_ungapped")
    if not col_ident_gaps:
        hidden.append("ident_with_gaps")
    if not col_indel:
        hidden.append("indel")
    if not col_region:
        hidden.append("region")
    return {
        "no_title": not fig_title,
        "no_subtitle": not fig_subtitle,
        "no_gap_track": not gap_track,
        "no_gene_track": not gene_track,
        "no_footer": not fig_footer,
        "no_table_title": not table_title,
        "no_table_subtitle": not table_subtitle,
        "no_table_note": not table_note,
        "hide_table_column": hidden,
    }


def parse_aliases(pairs: Optional[Sequence[str]]) -> Dict[str, str]:
    """Turn ``["old=new", ...]`` into a mapping, keeping any ``=`` inside the alias."""
    aliases: Dict[str, str] = {}
    for item in pairs or ():
        if "=" not in item:
            raise ValueError(f"--alias expects NAME=ALIAS, got {item!r}")
        old, new = item.split("=", 1)
        old, new = old.strip(), new.strip()
        if not old or not new:
            raise ValueError(f"--alias expects NAME=ALIAS, got {item!r}")
        aliases[old] = new
    return aliases


def sequence_aliases(opts) -> Dict[str, str]:
    aliases = getattr(opts, "aliases", None)
    if aliases is not None:
        return aliases
    return parse_aliases(getattr(opts, "alias", None))


class DuplicateDisplayNameError(ValueError):
    """Raised when aliases would make two sequences share a display name."""

    def __init__(self, names: Sequence[str]):
        self.names = list(names)
        super().__init__(
            "duplicate display name(s) after applying aliases: " + ", ".join(self.names)
        )


def duplicate_display_names(names: Sequence[str], aliases: Dict[str, str]) -> List[str]:
    """Display names that would appear more than once.

    A sequence keeps its original name unless an alias is given, so a clash can
    come from aliasing one sequence onto another's name or onto a shared alias.
    """
    counts: Dict[str, int] = {}
    for name in names:
        shown = aliases.get(name, name)
        counts[shown] = counts.get(shown, 0) + 1
    return sorted(shown for shown, count in counts.items() if count > 1)


def display_names(result: AnalysisResult, opts) -> List[str]:
    """Sequence names as they should appear in figures, tables and Markdown."""
    aliases = sequence_aliases(opts)
    return [aliases.get(name, name) for name in result.alignment.names]


def display_name(result: AnalysisResult, index: int, opts) -> str:
    aliases = sequence_aliases(opts)
    return aliases.get(result.alignment.names[index], result.alignment.names[index])


def unused_aliases(result: AnalysisResult, opts) -> List[str]:
    known = set(result.alignment.names)
    return [name for name in sequence_aliases(opts) if name not in known]


def hotspot_region(result: AnalysisResult, hotspot: Hotspot) -> str:
    """Name the annotated feature(s) overlapping a hotspot."""
    if not result.features:
        return "—"
    overlaps = []
    for feat in result.features:
        overlap = min(hotspot.end, feat.end) - max(hotspot.start, feat.start) + 1
        if overlap > 0:
            overlaps.append((overlap, feat.name))
    if not overlaps:
        return "—"
    overlaps.sort(key=lambda item: (-item[0], item[1]))
    best_overlap, best_name = overlaps[0]
    names = [best_name]
    if len(overlaps) > 1:
        second_overlap, second_name = overlaps[1]
        if second_name != best_name and second_overlap >= 0.2 * hotspot.length:
            names.append(second_name)
    return " / ".join(names)


def hotspot_rows(result: AnalysisResult, lang: str = "en") -> List[Dict[str, object]]:
    """Rows for the hotspot table shared by the GUI, Markdown and PDF output."""
    del lang  # formatting is language independent; headers are not
    rows: List[Dict[str, object]] = []
    for index, hs in enumerate(result.hotspots, start=1):
        rows.append(
            {
                "index": index,
                "start": hs.start,
                "end": hs.end,
                "coords": f"{hs.start:,}–{hs.end:,}",
                "length": hs.length,
                "ident_ungapped": round(hs.ident, 2),
                "ident_with_gaps": round(hs.ident_with_gaps, 2),
                "indel": f"{hs.gap_frac:.1f} / {hs.indel_bp_pct:.1f}",
                "indel_gap_pct": round(hs.gap_frac, 2),
                "indel_bp_pct": round(hs.indel_bp_pct, 2),
                "indel_bp": hs.indel_bp,
                "region": hotspot_region(result, hs),
            }
        )
    return rows


def window_rows(result: AnalysisResult) -> List[Dict[str, object]]:
    rows = []
    for w in result.windows:
        rows.append(
            {
                "start": w.start,
                "end": w.end,
                "columns": w.ncols,
                "matches": w.matches,
                "mismatches": w.mismatches,
                "gaps": w.gaps,
                "identity_ungapped": round(w.ident_ungapped, 3),
                "identity_with_gaps": round(w.ident_with_gaps, 3),
                "gap_pct": round(w.gap_frac, 3),
            }
        )
    return rows


def write_windows_tsv(result: AnalysisResult, path: str) -> None:
    rows = window_rows(result)
    cols = list(rows[0].keys()) if rows else [
        "start", "end", "columns", "matches", "mismatches", "gaps",
        "identity_ungapped", "identity_with_gaps", "gap_pct",
    ]
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(cols) + "\n")
        for r in rows:
            handle.write("\t".join(str(r[c]) for c in cols) + "\n")


def hotspot_context(result: AnalysisResult, opts, lang: str = "en") -> str:
    aln = result.alignment
    shown = display_names(result, opts)
    names = f"{shown[0]} vs {shown[1] if aln.nseq > 1 else ''}"
    span = result.coord_max - result.coord_min + 1
    with_genes = not getattr(opts, "no_gene_track", False)
    if lang == "zh":
        text = (
            f"{names} ｜ 窗口 {opts.window} bp / 步长 {opts.step} bp ｜ "
            f"比较区间 {result.coord_min:,}–{result.coord_max:,}（{span:,} bp）｜ "
            f"一致度 {result.identity_ungapped:.1f}%"
        )
        if with_genes:
            text += f" ｜ 基因来源 {feature_source_label(result, 'zh')}"
        return text
    text = (
        f"{names} | window {opts.window} bp / step {opts.step} bp | "
        f"compared {result.coord_min:,}-{result.coord_max:,} ({span:,} bp) | "
        f"identity {result.identity_ungapped:.1f}%"
    )
    if with_genes:
        text += f" | genes: {feature_source_label(result, 'en')}"
    return text


def table_cell(row: Dict[str, object], column: str) -> str:
    """Format one hotspot row cell for the text and figure tables."""
    if column == "index":
        return str(row["index"])
    if column == "coords":
        return str(row["coords"])
    if column == "length":
        return str(row["length"])
    if column == "ident_ungapped":
        return f"{row['ident_ungapped']:.2f}"
    if column == "ident_with_gaps":
        return f"{row['ident_with_gaps']:.2f}"
    if column == "indel":
        return str(row["indel"])
    if column == "region":
        return str(row["region"])
    raise KeyError(column)


def write_hotspot_markdown(result: AnalysisResult, path: str, opts, lang: str = "en") -> None:
    headers_all = HOTSPOT_TABLE_HEADERS.get(lang, HOTSPOT_TABLE_HEADERS["en"])
    header_of = dict(zip(HOTSPOT_TABLE_COLUMNS, headers_all))
    columns = visible_table_columns(opts)
    headers = [header_of[col] for col in columns]
    rows = hotspot_rows(result, lang)
    aln = result.alignment
    shown = display_names(result, opts)

    lines: List[str] = []
    if not getattr(opts, "no_table_title", False):
        title = "低保守区" if lang == "zh" else "Low-conservation regions"
        lines.extend([f"# {title}", ""])

    if not getattr(opts, "no_table_subtitle", False):
        lines.append(f"- " + ("比对" if lang == "zh" else "Alignment")
                     + f": `{shown[0]}` vs `{shown[1] if aln.nseq > 1 else ''}`")
        lines.append(
            "- " + ("坐标序列" if lang == "zh" else "Coordinate sequence")
            + f": `{shown[result.coords_seq]}`"
        )
        lines.append(
            "- " + ("比较区间" if lang == "zh" else "Compared region")
            + f": {result.coord_min:,}–{result.coord_max:,} "
            f"({result.coord_max - result.coord_min + 1:,} bp)"
        )
        lines.append(
            "- " + ("窗口 / 步长" if lang == "zh" else "Window / step")
            + f": {opts.window} bp / {opts.step} bp"
        )
        lines.append(
            "- " + ("总体一致度" if lang == "zh" else "Overall identity")
            + f": {result.identity_ungapped:.1f}% "
            + ("(不计 gap)" if lang == "zh" else "(gap-free)")
        )
        if not getattr(opts, "no_gene_track", False):
            lines.append(
                "- " + ("基因来源" if lang == "zh" else "Gene source")
                + f": {feature_source_label(result, lang)}"
            )
        lines.append("")

    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")
    if rows:
        for row in rows:
            lines.append("| " + " | ".join(table_cell(row, col) for col in columns) + " |")
    else:
        lines.append("| " + " | ".join(["—"] * len(headers)) + " |")

    if not getattr(opts, "no_table_note", False):
        lines.append("")
        lines.append(HOTSPOT_TABLE_NOTE.get(lang, HOTSPOT_TABLE_NOTE["en"]))

    lines.append("")
    lines.append(
        "```bash\n"
        f"python align_conservation_plot.py -i <alignment> -o <figure.svg> "
        f"--coords-seq {result.coords_seq + 1} --window {opts.window} --step {opts.step}\n"
        "```"
    )
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def build_hotspot_table_figures(
    result: AnalysisResult, opts, lang: str = "en", rows_per_page: int = 32
) -> List["object"]:
    """Render the hotspot table as one or more matplotlib pages."""
    from matplotlib.figure import Figure

    headers_all = HOTSPOT_TABLE_HEADERS.get(lang, HOTSPOT_TABLE_HEADERS["en"])
    header_of = dict(zip(HOTSPOT_TABLE_COLUMNS, headers_all))
    columns = visible_table_columns(opts)
    headers = [header_of[col] for col in columns]
    rows = hotspot_rows(result, lang)
    heading = "低保守区一览" if lang == "zh" else "Low-conservation regions"
    subtitle = hotspot_context(result, opts, lang)
    note = HOTSPOT_TABLE_NOTE.get(lang, HOTSPOT_TABLE_NOTE["en"])
    show_title = not getattr(opts, "no_table_title", False)
    show_subtitle = not getattr(opts, "no_table_subtitle", False)
    show_note = not getattr(opts, "no_table_note", False)

    # Column widths follow the same proportions, renormalised to the visible set.
    base_widths = {
        "index": 0.06, "coords": 0.15, "length": 0.10,
        "ident_ungapped": 0.17, "ident_with_gaps": 0.17, "indel": 0.20, "region": 0.15,
    }
    weights = [base_widths.get(col, 0.12) for col in columns]
    total_weight = sum(weights) or 1.0
    col_widths = [w / total_weight for w in weights]

    pages: List["object"] = []
    chunks = [rows[i : i + rows_per_page] for i in range(0, len(rows), rows_per_page)] or [[]]
    for page_index, chunk in enumerate(chunks, start=1):
        fig = Figure(figsize=(11.0, 6.6), dpi=opts.dpi, layout="constrained")
        ax = fig.add_subplot(111)
        ax.axis("off")

        if show_title:
            title = heading if len(chunks) == 1 else f"{heading} ({page_index}/{len(chunks)})"
            ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=14)
        if show_subtitle:
            ax.text(
                0.0, 1.0, subtitle, transform=ax.transAxes, fontsize=7.5, color="#555",
                va="bottom", ha="left",
            )

        table_rows = [list(headers)] + [
            [table_cell(r, col) for col in columns] for r in chunk
        ]
        if not chunk:
            table_rows.append(["—"] * len(headers))

        top = 0.90 if (show_title or show_subtitle) else 0.98
        n_rows = len(table_rows)
        row_h = min(0.075, (top - 0.10) / max(1, n_rows))
        table_h = row_h * n_rows
        table = ax.table(
            cellText=table_rows[1:],
            colLabels=table_rows[0],
            cellLoc="center",
            colLoc="center",
            bbox=[0.0, top - table_h, 1.0, table_h],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8.5 if n_rows > 12 else 9.5)
        for (row_idx, col_idx), cell in table.get_celld().items():
            cell.set_edgecolor("#cccccc")
            if col_idx < len(col_widths):
                cell.set_width(col_widths[col_idx])
            if row_idx == 0:
                cell.set_facecolor("#eef3f8")
                cell.set_text_props(fontweight="bold", fontsize=8.5)
            elif row_idx % 2 == 0:
                cell.set_facecolor("#fafafa")

        if show_note:
            ax.text(0.0, 0.02, note, transform=ax.transAxes, fontsize=7, color="#666", va="bottom")
        pages.append(fig)
    return pages


def save_pdf_bundle(
    fig: "object",
    result: AnalysisResult,
    opts,
    path: str,
    include_figure: bool = True,
    include_table: bool = True,
    rows_per_page: int = 32,
    lang: Optional[str] = None,
) -> int:
    """Write a PDF with the figure first and the hotspot table on later pages.

    Returns the number of pages written.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    pages = 0
    with PdfPages(path) as pdf:
        if include_figure and fig is not None:
            pdf.savefig(fig, bbox_inches="tight")
            pages += 1
        if include_table:
            for table_fig in build_hotspot_table_figures(
                result, opts, lang or getattr(opts, "lang", "en"), rows_per_page
            ):
                pdf.savefig(table_fig, bbox_inches="tight")
                pages += 1
    return pages


def summary_dict(result: AnalysisResult, opts) -> Dict[str, object]:
    aln = result.alignment
    aliases = sequence_aliases(opts)
    return {
        "sequences": aln.names,
        "sequence_aliases": aliases,
        "display_names": display_names(result, opts),
        "source_format": aln.source_format,
        "alignment_columns": aln.ncol,
        "coordinate_sequence": aln.names[result.coords_seq],
        "compared_region": [result.coord_min, result.coord_max],
        "compared_bp": result.coord_max - result.coord_min + 1,
        "terminal_all_gap_columns": [result.lead_gap, result.trail_gap],
        "core_columns": result.core_columns,
        "matches": result.total_matches,
        "mismatches": result.total_mismatches,
        "gap_columns": result.total_gaps,
        "identity_ungapped_pct": round(result.identity_ungapped, 3),
        "identity_with_gaps_pct": round(result.identity_with_gaps, 3),
        "window": opts.window,
        "step": opts.step,
        "indel_threshold_bp": opts.indel_min,
        "indels": [
            {
                "sequence": aln.names[i.seq_index],
                "kind": i.kind,
                "position": i.coord,
                "length": i.length,
            }
            for i in result.indels
        ],
        "hotspots": [
            {
                "index": row["index"],
                "start": row["start"],
                "end": row["end"],
                "length": row["length"],
                "identity_pct": row["ident_ungapped"],
                "identity_with_gaps_pct": row["ident_with_gaps"],
                "gap_pct": row["indel_gap_pct"],
                "indel_bp": row["indel_bp"],
                "indel_bp_pct": row["indel_bp_pct"],
                "region": row["region"],
                "columns": h.ncols,
            }
            for row, h in zip(hotspot_rows(result, getattr(opts, "lang", "en")), result.hotspots)
        ],
        "annotation_source": result.feature_source,
        "annotation_origin": getattr(result, "feature_source_origin", "") or None,
        "features": [
            {"name": f.name, "start": f.start, "end": f.end, "strand": f.strand}
            for f in result.features
        ],
        "warnings": result.warnings,
    }


def _fmt_num(value: float) -> str:
    return f"{value:,.0f}" if float(value).is_integer() else f"{value:,.1f}"


def print_report(result: AnalysisResult, opts) -> None:
    aln = result.alignment
    aliases = sequence_aliases(opts)
    shown = display_names(result, opts)
    sequence_label = ", ".join(
        f"{display} ({original})" if display != original else display
        for display, original in zip(shown, aln.names)
    )
    print("=" * 72)
    print(f"sequences        : {sequence_label}")
    print(f"format           : {aln.source_format}")
    print(f"alignment columns: {aln.ncol:,}")
    print(
        f"terminal gaps    : {result.lead_gap:,} (5') / {result.trail_gap:,} (3')"
    )
    print(
        f"compared         : {shown[result.coords_seq]} "
        f"{result.coord_min:,}-{result.coord_max:,} ({result.coord_max - result.coord_min + 1:,} bp)"
    )
    print(
        f"identity         : {result.identity_ungapped:.1f}% (gap-free) / "
        f"{result.identity_with_gaps:.1f}% (gaps counted)"
    )
    print(f"matches          : {_fmt_num(result.total_matches)}")
    print(f"mismatches       : {_fmt_num(result.total_mismatches)}")
    print(f"gap columns      : {result.total_gaps:,}")
    print(f"genome median    : {statistics_median([w.ident_ungapped for w in result.windows if not math.isnan(w.ident_ungapped)]):.1f}%")
    if result.hotspots:
        print("-" * 72)
        print("low-conservation hotspots:")
        for i, h in enumerate(result.hotspots[: opts.top], start=1):
            print(
                f"  {i:2d}. {h.start:>7,}-{h.end:<7,}  {h.end - h.start + 1:>5,} bp  "
                f"identity {h.ident:5.1f}% (with gaps {h.ident_with_gaps:5.1f}%)  "
                f"indel {h.gap_frac:4.1f}%/{h.indel_bp_pct:4.1f}%  "
                f"{hotspot_region(result, h)}"
            )
    conserved = sorted(
        [w for w in result.windows if not math.isnan(w.ident_ungapped) and w.ncols >= opts.window * 0.5],
        key=lambda w: -w.ident_ungapped,
    )[: opts.top]
    if conserved:
        print("-" * 72)
        print("most conserved windows:")
        for i, w in enumerate(conserved, start=1):
            print(
                f"  {i:2d}. {w.start:>7,}-{w.end:<7,}  identity {w.ident_ungapped:5.1f}%"
            )
    print("=" * 72)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Sliding-window conservation plot for sequence alignments.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", help="alignment file (srspair/pair/FASTA/Clustal/...)")
    p.add_argument("-o", "--out", help="output figure path (.svg)")
    p.add_argument("--pdf", nargs="?", const="", default=None, help="also write PDF (optional path)")
    p.add_argument("--png", nargs="?", const="", default=None, help="also write PNG (optional path)")
    p.add_argument(
        "--annotation",
        help="GenBank / GFF3 / TSV annotation file, or an NCBI accession "
        "(e.g. NC_004718.3) to download on the fly",
    )
    p.add_argument("--annotation-format", choices=["genbank", "gff3", "tsv"], help="override format detection")
    p.add_argument("--record", help="GenBank record id to use for annotation")
    p.add_argument(
        "--annotation-timeout",
        type=float,
        default=60.0,
        help="seconds to wait for an online annotation download",
    )
    p.add_argument(
        "--no-annotation-cache",
        action="store_true",
        help="do not read or write the local download cache",
    )
    p.add_argument(
        "--refresh-annotation",
        action="store_true",
        help="ignore the cached copy and download the accession again",
    )
    p.add_argument("--min-orf-aa", type=int, default=90, help="minimum ORF length when predicting genes")
    p.add_argument("--coords-seq", type=int, default=1, help="1-based index of the coordinate reference sequence")
    p.add_argument("--window", type=int, default=300, help="window size in bp")
    p.add_argument("--step", type=int, default=50, help="step size in bp")
    p.add_argument("--indel-min", type=int, default=8, help="minimum indel length to report")
    p.add_argument("--hotspot-k", type=float, default=1.5, help="hotspot threshold = median - k * 1.4826 * MAD")
    p.add_argument("--hotspot-pct", type=float, default=None, help="use this percentile instead of MAD rule")
    p.add_argument("--keep-terminal-gaps", action="store_true", help="keep terminal all-gap columns")
    p.add_argument("--title", help="figure title")
    p.add_argument("--lang", choices=["en", "zh"], default="en", help="figure label language")
    p.add_argument("--figsize", nargs=2, type=float, metavar=("W", "H"), help="figure size in inches")
    p.add_argument("--dpi", type=int, default=300, help="raster/export resolution")
    p.add_argument("--stats-out", help="write per-window statistics TSV")
    p.add_argument("--summary-out", help="write JSON summary")
    p.add_argument("--hotspot-table-md", help="write the low-conservation table as Markdown")
    p.add_argument(
        "--pdf-with-table",
        action="store_true",
        help="append the hotspot table pages to the PDF (figure stays on page 1)",
    )
    p.add_argument(
        "--alias",
        action="append",
        default=[],
        metavar="NAME=ALIAS",
        help="display alias for a sequence name (repeatable)",
    )
    p.add_argument("--no-title", action="store_true", help="hide the figure title")
    p.add_argument("--no-subtitle", action="store_true", help="hide the figure subtitle")
    p.add_argument("--no-gap-track", action="store_true", help="hide the gap%% track")
    p.add_argument(
        "--no-gene-track",
        action="store_true",
        help="hide the gene track; annotation is not loaded and the Region column is dropped",
    )
    p.add_argument("--no-footer", action="store_true", help="hide the figure footer note")
    p.add_argument("--no-table-title", action="store_true", help="hide the table title")
    p.add_argument("--no-table-subtitle", action="store_true", help="hide the table subtitle")
    p.add_argument(
        "--hide-table-column",
        action="append",
        default=[],
        choices=list(HIDEABLE_TABLE_COLUMNS),
        help="hide a table column (repeatable)",
    )
    p.add_argument(
        "--no-table-note",
        action="store_true",
        help="hide the table note and the Indel %% column together",
    )
    p.add_argument(
        "--show-font",
        action="store_true",
        help="print the figure font and its licence, then exit",
    )
    p.add_argument("--top", type=int, default=10, help="number of hotspots to print")
    p.add_argument("--gui", action="store_true", help="launch the Tkinter GUI")
    p.add_argument("--version", action="version", version="align_conservation_plot 1.0")
    return p


def render_and_save(result: AnalysisResult, opts, svg_path: Optional[str] = None) -> "object":
    fig = build_figure(result, opts)
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    FigureCanvasAgg(fig)
    if svg_path:
        fig.savefig(svg_path, format="svg", bbox_inches="tight")
    if opts.pdf is not None:
        pdf_path = opts.pdf or os.path.splitext(svg_path or "figure")[0] + ".pdf"
        pages = save_pdf_bundle(
            fig,
            result,
            opts,
            pdf_path,
            include_figure=True,
            include_table=bool(getattr(opts, "pdf_with_table", False)),
        )
        print(f"PDF written: {pdf_path} ({pages} page{'s' if pages != 1 else ''})")
    if opts.png is not None:
        png_path = opts.png or os.path.splitext(svg_path or "figure")[0] + ".png"
        fig.savefig(png_path, format="png", dpi=opts.dpi, bbox_inches="tight")
        print(f"PNG written: {png_path}")
    return fig


def run_cli(args) -> int:
    if not args.input:
        print("error: --input is required (or use --gui)", file=sys.stderr)
        return 2
    out = args.out or (os.path.splitext(args.input)[0] + "_conservation.svg")

    try:
        return _run_cli_inner(args, out)
    except FileNotFoundError as exc:
        print(f"error: file not found - {exc.filename}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        if os.environ.get("ALIGN_PLOT_DEBUG"):
            traceback.print_exc()
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _run_cli_inner(args, out: str) -> int:
    configure_matplotlib(args.lang)
    # Validate the alias syntax before doing any expensive work.
    args.aliases = parse_aliases(getattr(args, "alias", None))
    aln = read_alignment(args.input)
    clashes = duplicate_display_names(aln.names, args.aliases)
    if clashes:
        raise DuplicateDisplayNameError(clashes)
    result = analyse(
        aln,
        coords_seq=args.coords_seq,
        window=args.window,
        step=args.step,
        indel_min=args.indel_min,
        hotspot_k=args.hotspot_k,
        hotspot_pct=args.hotspot_pct,
        keep_terminal_gaps=args.keep_terminal_gaps,
    )
    if getattr(args, "no_gene_track", False):
        # The track is hidden, so do not read the annotation at all.
        result.features = []
        result.feature_source = "none"
    else:
        attach_features(
            result,
            args.annotation,
            args.annotation_format,
            args.record,
            args.min_orf_aa,
            annotation_options={
                "timeout": args.annotation_timeout,
                "use_cache": not args.no_annotation_cache,
                "refresh": args.refresh_annotation,
            },
        )
    for name in unused_aliases(result, args):
        print(f"warning: alias for unknown sequence {name!r} ignored", file=sys.stderr)
    for warning in result.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    print_report(result, args)
    render_and_save(result, args, out)
    print(f"SVG written: {out}")
    if args.stats_out:
        write_windows_tsv(result, args.stats_out)
        print(f"window statistics: {args.stats_out}")
    if args.summary_out:
        with open(args.summary_out, "w", encoding="utf-8") as handle:
            json.dump(summary_dict(result, args), handle, ensure_ascii=False, indent=2)
        print(f"summary: {args.summary_out}")
    if args.hotspot_table_md:
        write_hotspot_markdown(result, args.hotspot_table_md, args, args.lang)
        print(f"hotspot table: {args.hotspot_table_md}")
    return 0


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------


GUI_TEXT = {
    "zh": {
        "app_title": "比对保守度可视化",
        "sec_input": "输入",
        "lbl_input": "比对文件",
        "btn_browse": "浏览…",
        "lbl_annotation": "注释文件或 accession（可空）",
        "hint_annotation": "留空则自动预测 ORF；也可直接填 NCBI Accession No.（形如 NC_004718.3）联网下载（结果会缓存）",
        "sec_params": "参数",
        "lbl_coords": "坐标参考序列 (1-based)",
        "lbl_window": "窗口大小 (bp)",
        "lbl_step": "步长 (bp)",
        "lbl_indel_min": "插入或删除 阈值 (bp)",
        "lbl_hotspot_k": "热点 k (MAD)",
        "lbl_hotspot_pct": "热点百分位 (空=用 k)",
        "lbl_title": "标题（可空）",
        "chk_trim": "剔除两端全 gap 列",
        "chk_gene_track": "基因轨道",
        "sec_alias": "序列别名（每行 原名=别名，可留空）",
        "sec_display": "显示开关",
        "chk_fig_title": "图标题",
        "chk_fig_subtitle": "图副标题",
        "chk_gap_track": "gap 轨道",
        "chk_footer": "图注释行",
        "chk_table_title": "表标题",
        "chk_table_subtitle": "表副标题",
        "chk_table_note": "表注释行",
        "sec_table_cols": "表格列",
        "chk_col_length": "长度",
        "chk_col_ident_ungapped": "一致度(不计gap)",
        "chk_col_ident_gaps": "一致度(计gap)",
        "chk_col_indel": "插入或删除 %",
        "chk_col_region": "所在区域",
        "hint_alias": "选好比对文件后会自动列出序列名；把等号右边改成想要的显示名即可。",
        "msg_alias_dup_title": "别名冲突",
        "msg_alias_dup_body": "应用别名后有多个序列使用同一个显示名：{names}。请修改别名后重试。",
        "lbl_fig_lang": "图表文字：",
        "lang_zh": "中文",
        "lang_en": "English",
        "btn_switch_lang": "Switch to English",
        "sec_export": "导出内容",
        "chk_data": "数据（窗口统计）",
        "lbl_data_format": "格式",
        "chk_svg": "图 SVG",
        "chk_png": "图 PNG",
        "lbl_dpi": "DPI",
        "chk_md": "低保守区表 Markdown",
        "chk_pdf": "PDF",
        "chk_pdf_fig": "图",
        "chk_pdf_table": "表",
        "btn_export": "一键导出…",
        "btn_draw": "绘图",
        "tab_figure": "图",
        "tab_table": "表",
        "placeholder": "请选择比对文件后点击「绘图」",
        "status_ready": "就绪",
        "status_running": "正在解析并计算…",
        "status_done": "完成：{nseq} 条序列，{cols:,} 列，一致度 {ident:.1f}%，热点 {hotspots} 处",
        "status_failed": "失败：{msg}",
        "status_exported": "已导出 {count} 个文件到 {path}",
        "msg_missing_title": "缺少输入",
        "msg_missing_body": "请先选择一个存在的比对文件。",
        "msg_param_title": "参数错误",
        "msg_error_title": "绘图失败",
        "msg_export_failed_title": "导出失败",
        "msg_export_title": "导出",
        "msg_export_none": "请至少勾选一项导出内容。",
        "msg_export_no_result": "请先点击「绘图」生成结果。",
        "msg_dpi_invalid": "DPI 必须是 72–1200 之间的整数。",
        "msg_dir_title": "选择导出目录",
        "msg_export_done_title": "导出完成",
        "msg_export_done_body": "已将 {files} 导出至：\n{path}",
        "export_list_sep": "、",
        "dialog_pick_input": "选择比对文件",
        "dialog_pick_annotation": "选择注释文件",
        "dialog_files_all": "所有文件",
        "table_empty": "尚未生成结果",
        "tree_header_index": "序号",
        "tree_header_region": "所在区域",
    },
    "en": {
        "app_title": "Alignment Conservation Plot",
        "sec_input": "Input",
        "lbl_input": "Alignment file",
        "btn_browse": "Browse…",
        "lbl_annotation": "Annotation file or accession (optional)",
        "hint_annotation": "Leave empty to predict ORFs; you may also paste an NCBI accession (e.g. NC_004718.3) to download it (cached).",
        "sec_params": "Parameters",
        "lbl_coords": "Coordinate sequence (1-based)",
        "lbl_window": "Window size (bp)",
        "lbl_step": "Step (bp)",
        "lbl_indel_min": "Indel threshold (bp)",
        "lbl_hotspot_k": "Hotspot k (MAD)",
        "lbl_hotspot_pct": "Hotspot percentile (empty = k)",
        "lbl_title": "Title (optional)",
        "chk_trim": "Trim terminal all-gap columns",
        "chk_gene_track": "Gene track",
        "sec_alias": "Sequence aliases (one name=alias per line, blanks ignored)",
        "sec_display": "Display switches",
        "chk_fig_title": "Figure title",
        "chk_fig_subtitle": "Figure subtitle",
        "chk_gap_track": "Gap track",
        "chk_footer": "Figure note",
        "chk_table_title": "Table title",
        "chk_table_subtitle": "Table subtitle",
        "chk_table_note": "Table note",
        "sec_table_cols": "Table columns",
        "chk_col_length": "Length",
        "chk_col_ident_ungapped": "Identity (ungapped)",
        "chk_col_ident_gaps": "Identity (with gaps)",
        "chk_col_indel": "Indel %",
        "chk_col_region": "Region",
        "hint_alias": "Sequence names appear once a file is chosen; edit the right-hand side to rename them.",
        "msg_alias_dup_title": "Alias conflict",
        "msg_alias_dup_body": "After applying the aliases, these display names are shared by more than one sequence: {names}. Change them and try again.",
        "lbl_fig_lang": "Figure text:",
        "lang_zh": "中文",
        "lang_en": "English",
        "btn_switch_lang": "切换到中文",
        "sec_export": "Export",
        "chk_data": "Data (window statistics)",
        "lbl_data_format": "Format",
        "chk_svg": "Figure SVG",
        "chk_png": "Figure PNG",
        "lbl_dpi": "DPI",
        "chk_md": "Hotspot table Markdown",
        "chk_pdf": "PDF",
        "chk_pdf_fig": "figure",
        "chk_pdf_table": "table",
        "btn_export": "Export selected…",
        "btn_draw": "Plot",
        "tab_figure": "Figure",
        "tab_table": "Table",
        "placeholder": "Choose an alignment file, then click Plot.",
        "status_ready": "Ready",
        "status_running": "Parsing and computing…",
        "status_done": "Done: {nseq} sequences, {cols:,} columns, identity {ident:.1f}%, {hotspots} hotspots",
        "status_failed": "Failed: {msg}",
        "status_exported": "Exported {count} file(s) to {path}",
        "msg_missing_title": "Missing input",
        "msg_missing_body": "Please choose an existing alignment file first.",
        "msg_param_title": "Invalid parameters",
        "msg_error_title": "Plot failed",
        "msg_export_failed_title": "Export failed",
        "msg_export_title": "Export",
        "msg_export_none": "Select at least one export option.",
        "msg_export_no_result": "Click Plot before exporting.",
        "msg_dpi_invalid": "DPI must be an integer between 72 and 1200.",
        "msg_dir_title": "Choose export folder",
        "msg_export_done_title": "Export complete",
        "msg_export_done_body": "Exported {files} to:\n{path}",
        "export_list_sep": ", ",
        "dialog_pick_input": "Choose alignment file",
        "dialog_pick_annotation": "Choose annotation file",
        "dialog_files_all": "All files",
        "table_empty": "No result yet",
        "tree_header_index": "#",
        "tree_header_region": "Region",
    },
}


def launch_gui(self_test: Optional[callable] = None) -> int:
    import queue
    import threading
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    from tkinter import font as tkfont

    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    configure_matplotlib("zh")

    root = tk.Tk()
    root.geometry("1380x880")
    root.minsize(1080, 680)
    ui_font = pick_ui_font(tkfont.families(root))
    try:
        tkfont.nametofont("TkDefaultFont").configure(family=ui_font, size=10)
        tkfont.nametofont("TkTextFont").configure(family=ui_font, size=10)
        tkfont.nametofont("TkMenuFont").configure(family=ui_font, size=10)
    except Exception:
        pass

    state = {
        "input": tk.StringVar(),
        "annotation": tk.StringVar(),
        "coords_seq": tk.StringVar(value="1"),
        "window": tk.StringVar(value="300"),
        "step": tk.StringVar(value="50"),
        "indel_min": tk.StringVar(value="8"),
        "hotspot_k": tk.StringVar(value="1.5"),
        "hotspot_pct": tk.StringVar(value=""),
        "title": tk.StringVar(),
        "lang": tk.StringVar(value="en"),
        "trim": tk.BooleanVar(value=True),
        "exp_data": tk.BooleanVar(value=True),
        "exp_data_format": tk.StringVar(value="TSV"),
        "exp_svg": tk.BooleanVar(value=True),
        "exp_png": tk.BooleanVar(value=False),
        "exp_md": tk.BooleanVar(value=True),
        "exp_pdf_fig": tk.BooleanVar(value=False),
        "exp_pdf_table": tk.BooleanVar(value=False),
        "png_dpi": tk.StringVar(value="300"),
        "gene_track": tk.BooleanVar(value=True),
        "fig_title": tk.BooleanVar(value=True),
        "fig_subtitle": tk.BooleanVar(value=True),
        "gap_track": tk.BooleanVar(value=True),
        "fig_footer": tk.BooleanVar(value=True),
        "table_title": tk.BooleanVar(value=True),
        "table_subtitle": tk.BooleanVar(value=True),
        "table_note": tk.BooleanVar(value=True),
        "col_length": tk.BooleanVar(value=True),
        "col_ident_ungapped": tk.BooleanVar(value=True),
        "col_ident_gaps": tk.BooleanVar(value=True),
        "col_indel": tk.BooleanVar(value=True),
        "col_region": tk.BooleanVar(value=True),
    }
    alias_text: "tk.Text"

    gui_lang = {"value": "en"}
    task_queue: "queue.Queue[tuple]" = queue.Queue()
    registry: List[Tuple[object, str, str]] = []
    app = {"fig": None, "result": None, "args": None, "canvas": None, "toolbar": None}

    def refresh_font_notice() -> None:
        """State which font renders the figures, and under which licence."""
        try:
            chosen = pick_figure_font(available_font_families())
        except Exception:  # noqa: BLE001 - a missing font must not break start-up
            chosen = None
        font_var.set(font_notice(chosen, gui_lang["value"]))

    def text(key: str) -> str:
        return GUI_TEXT[gui_lang["value"]].get(key, key)

    def register(widget, key: str, attr: str = "text"):
        registry.append((widget, key, attr))
        try:
            widget.configure(**{attr: text(key)})
        except tk.TclError:
            pass

    # The control panel scrolls: with every switch on it is taller than a laptop screen.
    left_outer = ttk.Frame(root, padding=(8, 8, 4, 8))
    left_outer.pack(side="left", fill="y")
    left_canvas = tk.Canvas(left_outer, width=478, highlightthickness=0, borderwidth=0)
    left_scroll = ttk.Scrollbar(left_outer, orient="vertical", command=left_canvas.yview)
    left_canvas.configure(yscrollcommand=left_scroll.set)
    left_scroll.pack(side="right", fill="y")
    left_canvas.pack(side="left", fill="y", expand=True)
    left = ttk.Frame(left_canvas, padding=(10, 4, 12, 8))
    left_window = left_canvas.create_window((0, 0), window=left, anchor="nw")
    def _fit_panel(_event=None):
        """Keep the scrollregion in sync; the width is fixed after measuring."""
        left_canvas.configure(scrollregion=left_canvas.bbox("all"))

    left.bind("<Configure>", _fit_panel)
    left_canvas.bind(
        "<Configure>",
        lambda e: left_canvas.itemconfigure(left_window, width=e.width),
    )

    def _scroll_pixels(pixels: float) -> None:
        """Scroll the panel by pixels; canvas "units" are too coarse for trackpads."""
        total = max(1, left.winfo_height())
        first, _last = left_canvas.yview()
        target = first * total + pixels
        left_canvas.yview_moveto(max(0.0, min(1.0, target / total)))

    def _in_panel(widget) -> bool:
        """True when widget sits inside the scrolling panel."""
        while widget is not None:
            if widget is alias_text:
                # Let the text box scroll itself.
                return False
            if widget is left_outer:
                return True
            parent = widget.winfo_parent()
            if not parent:
                return False
            try:
                widget = widget.nametowidget(parent)
            except KeyError:
                return False
        return False

    def _wheel_target(event):
        """Panel membership from the event itself, else from the pointer."""
        widget = getattr(event, "widget", None)
        if widget is not None:
            return _in_panel(widget)
        try:
            widget = root.winfo_containing(root.winfo_pointerx(), root.winfo_pointery())
        except tk.TclError:
            return False
        return widget is not None and _in_panel(widget)

    def _on_wheel(event):
        widget = getattr(event, "widget", None)
        trace = os.environ.get("ALIGN_PLOT_TRACE_WHEEL")
        if trace:
            print(
                f"wheel widget={widget} delta={getattr(event, 'delta', None)} "
                f"num={getattr(event, 'num', None)} in_panel={_in_panel(widget)} target={_wheel_target(event)}",
                file=sys.stderr,
            )
        if not _wheel_target(event):
            return None
        delta = getattr(event, "delta", 0)
        num = getattr(event, "num", 0)
        if delta:
            # macOS reports ~1 per notch (trackpads send fractions); Windows 120.
            pixels = -delta * 15.0 if abs(delta) <= 10 else -delta / 120.0 * 36.0
        elif num in (4, 5):
            pixels = -36.0 if num == 4 else 36.0
        else:
            return None
        _scroll_pixels(pixels)
        return "break"

    def _on_touchpad(event):
        """Tk 9 sends trackpad gestures as <TouchpadScroll>, not <MouseWheel>.

        %D packs the horizontal delta in the high 16 bits and the vertical
        delta in the low 16 bits, both signed (see tk::PreciseScrollDeltas).
        """
        widget = getattr(event, "widget", None)
        if os.environ.get("ALIGN_PLOT_TRACE_WHEEL"):
            print(
                f"touchpad widget={widget} delta={getattr(event, 'delta', None)} "
                f"in_panel={_in_panel(widget)} target={_wheel_target(event)}",
                file=sys.stderr,
            )
        if not _wheel_target(event):
            return None
        try:
            packed = int(getattr(event, "delta", 0)) & 0xFFFFFFFF
        except (TypeError, ValueError):
            return None
        dy = packed & 0xFFFF
        if dy >= 0x8000:
            dy -= 0x10000
        if not dy:
            return None
        _scroll_pixels(-float(dy))
        return "break"

    # A single binding on the toplevel is far more reliable than per-widget
    # bindings: every child inherits it, including widgets created later.
    root.bind_all("<MouseWheel>", _on_wheel)
    root.bind_all("<Button-4>", _on_wheel)
    root.bind_all("<Button-5>", _on_wheel)
    try:
        root.bind_all("<TouchpadScroll>", _on_touchpad)
        touchpad_events = True
    except tk.TclError:  # Tk < 9.1 has no trackpad event, wheels still work
        touchpad_events = False

    left.grid_columnconfigure(0, minsize=175)
    right = ttk.Frame(root, padding=(4, 8, 8, 8))
    right.pack(side="right", fill="both", expand=True)

    row = 0
    header = ttk.Frame(left)
    header.grid(row=row, column=0, columnspan=3, sticky="we", pady=(0, 2))
    title_label = ttk.Label(header, font=(ui_font, 11, "bold"))
    title_label.pack(side="left")
    register(title_label, "app_title")
    switch_btn = ttk.Button(header, width=16)
    switch_btn.pack(side="right")
    register(switch_btn, "btn_switch_lang")
    row += 1

    sec_input = ttk.Label(left, font=(ui_font, 9, "bold"))
    sec_input.grid(row=row, column=0, columnspan=3, sticky="w", pady=(4, 2))
    register(sec_input, "sec_input")
    row += 1

    lbl_input = ttk.Label(left, anchor="w")
    lbl_input.grid(row=row, column=0, columnspan=3, sticky="w")
    register(lbl_input, "lbl_input")
    row += 1
    ttk.Entry(left, textvariable=state["input"], width=26).grid(
        row=row, column=0, columnspan=2, sticky="we", pady=3
    )

    def pick_input():
        path = filedialog.askopenfilename(
            title=text("dialog_pick_input"),
            filetypes=[
                ("Alignment", "*.srspair *.pair *.aln *.fasta *.fa *.fas *.sto *.phy *.msf *.nex *.gz"),
                (text("dialog_files_all"), "*.*"),
            ],
        )
        if path:
            state["input"].set(path)
            fill_aliases()

    browse_in = ttk.Button(left, width=8, command=pick_input)
    browse_in.grid(row=row, column=2, padx=(4, 0))
    register(browse_in, "btn_browse")
    row += 1

    # The gene-track switch lives with the annotation input it controls.
    ann_head = ttk.Frame(left)
    ann_head.grid(row=row, column=0, columnspan=3, sticky="w", pady=(8, 0))
    gene_chk = ttk.Checkbutton(ann_head, variable=state["gene_track"])
    gene_chk.pack(side="left")
    register(gene_chk, "chk_gene_track")
    lbl_annotation = ttk.Label(ann_head)
    lbl_annotation.pack(side="left", padx=(8, 0))
    register(lbl_annotation, "lbl_annotation")
    row += 1
    anno_entry = ttk.Entry(left, textvariable=state["annotation"], width=26)
    anno_entry.grid(row=row, column=0, columnspan=2, sticky="we", pady=3)

    def pick_annotation():
        path = filedialog.askopenfilename(
            title=text("dialog_pick_annotation"),
            filetypes=[
                ("Annotation", "*.gb *.gbk *.genbank *.gff *.gff3 *.tsv *.txt"),
                (text("dialog_files_all"), "*.*"),
            ],
        )
        if path:
            state["annotation"].set(path)

    browse_ann = ttk.Button(left, width=8, command=pick_annotation)
    browse_ann.grid(row=row, column=2, padx=(4, 0))
    register(browse_ann, "btn_browse")
    row += 1
    hint = ttk.Label(left, foreground="#777", wraplength=320, font=(ui_font, 8))
    hint.grid(row=row, column=0, columnspan=3, sticky="w")
    register(hint, "hint_annotation")
    row += 1

    # --- sequence aliases --------------------------------------------------
    # Collapsed by default: most runs do not rename sequences.
    state["alias_open"] = tk.BooleanVar(value=False)
    alias_toggle = ttk.Checkbutton(left, variable=state["alias_open"], style="Toolbutton")
    alias_toggle.grid(row=row, column=0, columnspan=3, sticky="w", pady=(8, 0))
    register(alias_toggle, "sec_alias")
    row += 1
    alias_box = ttk.Frame(left)
    alias_box.grid(row=row, column=0, columnspan=3, sticky="we")
    alias_text = tk.Text(
        alias_box,
        height=5,
        width=20,
        wrap="none",
        font=(ui_font, 9),
        relief="solid",
        borderwidth=1,
        highlightthickness=0,
        padx=4,
        pady=2,
    )
    alias_scroll = ttk.Scrollbar(alias_box, orient="vertical", command=alias_text.yview)
    alias_text.configure(yscrollcommand=alias_scroll.set)
    alias_text.pack(side="left", fill="both", expand=True)
    alias_scroll.pack(side="right", fill="y")
    row += 1
    alias_hint = ttk.Label(left, foreground="#777", wraplength=320, font=(ui_font, 8))
    alias_hint.grid(row=row, column=0, columnspan=3, sticky="w")
    register(alias_hint, "hint_alias")
    row += 1

    def refresh_alias_visibility(*_args):
        """Collapse or expand the alias editor and keep the arrow in step."""
        open_ = state["alias_open"].get()
        alias_toggle.configure(text=("▾ " if open_ else "▸ ") + text("sec_alias"))
        if open_:
            alias_box.grid()
            alias_hint.grid()
        else:
            alias_box.grid_remove()
            alias_hint.grid_remove()
        root.after_idle(_fit_panel)

    state["alias_open"].trace_add("write", refresh_alias_visibility)

    def fill_aliases():
        path = state["input"].get().strip()
        names = read_sequence_names(path) if path else []
        alias_text.delete("1.0", "end")
        alias_text.insert("1.0", "\n".join(f"{name}=" for name in names))

    def read_alias_text() -> dict:
        pairs = []
        for line in alias_text.get("1.0", "end").splitlines():
            line = line.strip()
            if line and "=" in line:
                pairs.append(line)
        return parse_aliases(pairs)

    ttk.Separator(left, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="we", pady=4)
    row += 1
    sec_params = ttk.Label(left, font=(ui_font, 9, "bold"))
    sec_params.grid(row=row, column=0, columnspan=3, sticky="w")
    register(sec_params, "sec_params")
    row += 1

    def labelled_entry(parent, row_index, key, var, width=6, col=0, span=1):
        """One label/entry pair; two pairs fit side by side now the fonts are smaller."""
        label = ttk.Label(parent, anchor="w")
        label.grid(row=row_index, column=col, sticky="w", pady=1, padx=(0, 6))
        register(label, key)
        ttk.Entry(parent, textvariable=var, width=width).grid(
            row=row_index, column=col + 1, columnspan=span, sticky="we", pady=1, padx=(0, 12)
        )

    params_top = ttk.Frame(left)
    params_top.grid(row=row, column=0, columnspan=3, sticky="we")
    # Longest labels go in the first column of each row so the shared column
    # width stays small and the two pairs fit on one line.
    for pair_row, (left_key, left_var, right_key, right_var) in enumerate((
        ("lbl_coords", state["coords_seq"], "lbl_window", state["window"]),
        ("lbl_indel_min", state["indel_min"], "lbl_step", state["step"]),
        ("lbl_hotspot_pct", state["hotspot_pct"], "lbl_hotspot_k", state["hotspot_k"]),
    )):
        labelled_entry(params_top, pair_row, left_key, left_var, col=0)
        labelled_entry(params_top, pair_row, right_key, right_var, col=2)
    row += 1
    labelled_entry(left, row, "lbl_title", state["title"], width=20, span=2)
    row += 1

    trim_row = ttk.Frame(left)
    trim_row.grid(row=row, column=0, columnspan=3, sticky="w", pady=1)
    trim_chk = ttk.Checkbutton(trim_row, variable=state["trim"])
    trim_chk.pack(side="left")
    register(trim_chk, "chk_trim")
    fig_lang_label = ttk.Label(trim_row)
    fig_lang_label.pack(side="left", padx=(12, 0))
    register(fig_lang_label, "lbl_fig_lang")
    for value, key in (("zh", "lang_zh"), ("en", "lang_en")):
        rb = ttk.Radiobutton(trim_row, value=value, variable=state["lang"])
        rb.pack(side="left")
        register(rb, key)
    row += 1

    ttk.Separator(left, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="we", pady=4)
    row += 1

    # --- display switches --------------------------------------------------
    sec_display = ttk.Label(left, font=(ui_font, 9, "bold"))
    sec_display.grid(row=row, column=0, columnspan=3, sticky="w")
    register(sec_display, "sec_display")
    row += 1
    display_grid = ttk.Frame(left)
    display_grid.grid(row=row, column=0, columnspan=3, sticky="w")
    display_keys = (
        ("fig_title", "chk_fig_title"),
        ("fig_subtitle", "chk_fig_subtitle"),
        ("gap_track", "chk_gap_track"),
        ("fig_footer", "chk_footer"),
        ("table_title", "chk_table_title"),
        ("table_subtitle", "chk_table_subtitle"),
        ("table_note", "chk_table_note"),
    )
    display_checks = {}
    for i, (var_key, label_key) in enumerate(display_keys):
        chk = ttk.Checkbutton(display_grid, variable=state[var_key])
        chk.grid(row=i // 3, column=i % 3, sticky="w", padx=(0, 10))
        register(chk, label_key)
        display_checks[var_key] = chk
    row += 1

    # --- table columns -----------------------------------------------------
    sec_cols = ttk.Label(left, font=(ui_font, 9, "bold"))
    sec_cols.grid(row=row, column=0, columnspan=3, sticky="w", pady=(4, 0))
    register(sec_cols, "sec_table_cols")
    row += 1
    cols_grid = ttk.Frame(left)
    cols_grid.grid(row=row, column=0, columnspan=3, sticky="w")
    column_checks = {}
    for i, (var_key, label_key) in enumerate((
        ("col_length", "chk_col_length"),
        ("col_ident_ungapped", "chk_col_ident_ungapped"),
        ("col_ident_gaps", "chk_col_ident_gaps"),
        ("col_indel", "chk_col_indel"),
        ("col_region", "chk_col_region"),
    )):
        chk = ttk.Checkbutton(cols_grid, variable=state[var_key])
        chk.grid(row=i // 3, column=i % 3, sticky="w", padx=(0, 10))
        register(chk, label_key)
        column_checks[var_key] = chk
    row += 1

    region_last = {"value": True}

    def sync_gene_dependent_widgets():
        """Gene track off: no annotation input and no Region column."""
        on = state["gene_track"].get()
        anno_entry.configure(state="normal" if on else "disabled")
        browse_ann.configure(state="normal" if on else "disabled")
        region_chk = column_checks["col_region"]
        if on:
            region_chk.configure(state="normal")
            state["col_region"].set(region_last["value"])
        else:
            region_last["value"] = state["col_region"].get()
            state["col_region"].set(False)
            region_chk.configure(state="disabled")

    def sync_note_and_indel(*_args):
        """The table note and the Indel % column are shown or hidden together."""
        on = state["table_note"].get()
        state["col_indel"].set(on)
        column_checks["col_indel"].configure(state="normal" if on else "disabled")

    state["gene_track"].trace_add("write", lambda *_: sync_gene_dependent_widgets())
    state["table_note"].trace_add("write", sync_note_and_indel)
    def sync_note_from_indel(*_args):
        if state["col_indel"].get():
            state["table_note"].set(True)

    state["col_indel"].trace_add("write", sync_note_from_indel)

    ttk.Separator(left, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="we", pady=4)
    row += 1
    sec_export = ttk.Label(left, font=(ui_font, 9, "bold"))
    sec_export.grid(row=row, column=0, columnspan=3, sticky="w")
    register(sec_export, "sec_export")
    row += 1

    data_row = ttk.Frame(left)
    data_row.grid(row=row, column=0, columnspan=3, sticky="w")
    data_chk = ttk.Checkbutton(data_row, variable=state["exp_data"])
    data_chk.pack(side="left")
    register(data_chk, "chk_data")
    fmt_label = ttk.Label(data_row)
    fmt_label.pack(side="left", padx=(8, 2))
    register(fmt_label, "lbl_data_format")
    fmt_box = ttk.Combobox(
        data_row, textvariable=state["exp_data_format"], values=("TSV", "JSON"), width=6, state="readonly"
    )
    fmt_box.pack(side="left")
    row += 1

    svg_row = ttk.Frame(left)
    svg_row.grid(row=row, column=0, columnspan=3, sticky="w")
    svg_chk = ttk.Checkbutton(svg_row, variable=state["exp_svg"])
    svg_chk.pack(side="left")
    register(svg_chk, "chk_svg")
    png_chk = ttk.Checkbutton(svg_row, variable=state["exp_png"])
    png_chk.pack(side="left", padx=(12, 0))
    register(png_chk, "chk_png")
    dpi_label = ttk.Label(svg_row)
    dpi_label.pack(side="left", padx=(8, 2))
    register(dpi_label, "lbl_dpi")
    ttk.Entry(svg_row, textvariable=state["png_dpi"], width=5).pack(side="left")
    row += 1

    md_row = ttk.Frame(left)
    md_row.grid(row=row, column=0, columnspan=3, sticky="w")
    md_chk = ttk.Checkbutton(md_row, variable=state["exp_md"])
    md_chk.pack(side="left")
    register(md_chk, "chk_md")
    pdf_label = ttk.Label(md_row)
    pdf_label.pack(side="left", padx=(12, 0))
    register(pdf_label, "chk_pdf")
    for var, key in ((state["exp_pdf_fig"], "chk_pdf_fig"), (state["exp_pdf_table"], "chk_pdf_table")):
        sub = ttk.Checkbutton(md_row, variable=var)
        sub.pack(side="left", padx=(8, 0))
        register(sub, key)
    row += 1

    ttk.Separator(left, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="we", pady=4)
    row += 1

    status = tk.StringVar(value=text("status_ready"))
    action_row = ttk.Frame(left)
    action_row.grid(row=row, column=0, columnspan=3, sticky="we", pady=(0, 4))
    draw_btn = ttk.Button(action_row, width=10)
    draw_btn.pack(side="left")
    register(draw_btn, "btn_draw")
    export_btn = ttk.Button(action_row, width=14, state="disabled")
    export_btn.pack(side="left", padx=6)
    register(export_btn, "btn_export")
    row += 1
    font_var = tk.StringVar()
    font_label = ttk.Label(
        left, textvariable=font_var, foreground="#777", wraplength=460,
        font=(ui_font, 8), justify="left",
    )
    font_label.grid(row=row, column=0, columnspan=3, sticky="w", pady=(2, 2))
    row += 1
    status_label = ttk.Label(left, textvariable=status, foreground="#555", wraplength=340)
    status_label.grid(row=row, column=0, columnspan=3, sticky="w")

    notebook = ttk.Notebook(right)
    notebook.pack(fill="both", expand=True)
    tab_figure = ttk.Frame(notebook)
    tab_table = ttk.Frame(notebook)
    notebook.add(tab_figure, text=text("tab_figure"))
    notebook.add(tab_table, text=text("tab_table"))

    placeholder = ttk.Label(tab_figure, foreground="#888")
    placeholder.pack(expand=True)
    register(placeholder, "placeholder")

    table_head = ttk.Frame(tab_table)
    table_head.pack(fill="x", padx=8, pady=(8, 0))
    table_title_label = ttk.Label(table_head, font=(ui_font, 11, "bold"), anchor="w")
    table_title_label.pack(fill="x")
    table_subtitle_label = ttk.Label(table_head, foreground="#555", anchor="w",
                                     font=(ui_font, 8), wraplength=760)
    table_subtitle_label.pack(fill="x")

    table_holder = ttk.Frame(tab_table)
    table_holder.pack(fill="both", expand=True, padx=6, pady=6)
    tree = ttk.Treeview(
        table_holder, columns=HOTSPOT_TABLE_COLUMNS, show="headings", selectmode="browse"
    )
    tree_scroll = ttk.Scrollbar(table_holder, orient="vertical", command=tree.yview)
    tree.configure(yscrollcommand=tree_scroll.set)
    tree.pack(side="left", fill="both", expand=True)
    tree_scroll.pack(side="right", fill="y")
    widths = {"index": 45, "coords": 130, "length": 80, "ident_ungapped": 115,
              "ident_with_gaps": 115, "indel": 175, "region": 120}
    for col in HOTSPOT_TABLE_COLUMNS:
        tree.column(col, width=widths.get(col, 120), anchor="center", stretch=True)
    table_note = ttk.Label(tab_table, foreground="#666", wraplength=760, font=(ui_font, 9))
    table_note.pack(fill="x", padx=8, pady=(0, 8))
    empty_label = ttk.Label(tab_table, foreground="#888")
    register(empty_label, "table_empty")

    def table_lang() -> str:
        """The table belongs to the output, so it follows the figure language.

        The interface language only drives GUI chrome; otherwise setting the
        figure text to English while the interface stayed Chinese produced an
        English figure with a Chinese table.
        """
        return state["lang"].get() or "en"

    def refresh_table_headers():
        headers = HOTSPOT_TABLE_HEADERS.get(table_lang(), HOTSPOT_TABLE_HEADERS["en"])
        for col, header in zip(HOTSPOT_TABLE_COLUMNS, headers):
            tree.heading(col, text=header)
        table_note.configure(text=HOTSPOT_TABLE_NOTE.get(table_lang(), HOTSPOT_TABLE_NOTE["en"]))
        table_title_label.configure(
            text="低保守区一览" if table_lang() == "zh" else "Low-conservation regions"
        )

    def collect_display_opts() -> dict:
        """Display switches as parser-compatible option values."""
        return display_options(
            fig_title=state["fig_title"].get(),
            fig_subtitle=state["fig_subtitle"].get(),
            gap_track=state["gap_track"].get(),
            gene_track=state["gene_track"].get(),
            fig_footer=state["fig_footer"].get(),
            table_title=state["table_title"].get(),
            table_subtitle=state["table_subtitle"].get(),
            table_note=state["table_note"].get(),
            col_length=state["col_length"].get(),
            col_ident_ungapped=state["col_ident_ungapped"].get(),
            col_ident_gaps=state["col_ident_gaps"].get(),
            col_indel=state["col_indel"].get(),
            col_region=state["col_region"].get(),
        )

    def preview_opts():
        """Current display switches plus the analysis options behind the last draw."""
        base = app["args"]
        if base is None:
            try:
                base = current_args()
            except Exception:  # noqa: BLE001 - preview must never crash
                base = build_parser().parse_args([])
        opts_ns = argparse.Namespace(**vars(base))
        for key, value in collect_display_opts().items():
            setattr(opts_ns, key, value)
        return opts_ns

    def refresh_table_layout():
        """Apply column visibility, table title/subtitle and the note line."""
        result = app["result"]
        opts_like = preview_opts()
        visible = visible_table_columns(opts_like)
        tree.configure(displaycolumns=visible)
        show_title = not opts_like.no_table_title
        show_subtitle = not opts_like.no_table_subtitle and result is not None
        # Re-pack in a fixed order; "before=" refuses an unpacked reference.
        table_title_label.pack_forget()
        table_subtitle_label.pack_forget()
        if show_title:
            table_title_label.pack(fill="x")
        if show_subtitle:
            table_subtitle_label.configure(text=hotspot_context(result, opts_like, table_lang()))
            table_subtitle_label.pack(fill="x")
        if opts_like.no_table_note:
            table_note.pack_forget()
        else:
            table_note.pack(fill="x", padx=8, pady=(0, 8))

    def refresh_table_rows():
        tree.delete(*tree.get_children())
        result = app["result"]
        if result is None:
            empty_label.pack(expand=True)
            return
        empty_label.pack_forget()
        for row_data in hotspot_rows(result, table_lang()):
            tree.insert(
                "",
                "end",
                values=(
                    row_data["index"],
                    row_data["coords"],
                    f"{row_data['length']:,}",
                    f"{row_data['ident_ungapped']:.2f}",
                    f"{row_data['ident_with_gaps']:.2f}",
                    row_data["indel"],
                    row_data["region"],
                ),
            )

    def apply_language():
        for widget, key, attr in registry:
            try:
                widget.configure(**{attr: text(key)})
            except tk.TclError:
                pass
        root.title(text("app_title"))
        notebook.tab(0, text=text("tab_figure"))
        notebook.tab(1, text=text("tab_table"))
        switch_btn.configure(text=text("btn_switch_lang"))
        refresh_font_notice()
        refresh_alias_visibility()
        root.after_idle(_fit_panel)
        refresh_table_headers()
        refresh_table_rows()
        refresh_table_layout()
        if app["result"] is None:
            status.set(text("status_ready"))
        else:
            status.set(
                text("status_done").format(
                    nseq=app["result"].alignment.nseq,
                    cols=app["result"].core_columns,
                    ident=app["result"].identity_ungapped,
                    hotspots=len(app["result"].hotspots),
                )
            )

    def toggle_language():
        gui_lang["value"] = "en" if gui_lang["value"] == "zh" else "zh"
        apply_language()

    switch_btn.configure(command=toggle_language)
    state["lang"].trace_add("write", lambda *_args: refresh_table_headers())
    state["lang"].trace_add("write", lambda *_args: refresh_table_rows())
    state["lang"].trace_add("write", lambda *_args: refresh_table_layout())
    for _key in ("col_length", "col_ident_ungapped", "col_ident_gaps", "col_indel",
                 "col_region", "table_title", "table_subtitle", "table_note"):
        state[_key].trace_add("write", lambda *_args: refresh_table_layout())

    def current_args():
        args = build_parser().parse_args([])
        args.input = state["input"].get().strip()
        args.annotation = state["annotation"].get().strip() or None
        args.coords_seq = int(state["coords_seq"].get())
        args.window = int(state["window"].get())
        args.step = int(state["step"].get())
        args.indel_min = int(state["indel_min"].get())
        args.hotspot_k = float(state["hotspot_k"].get())
        pct = state["hotspot_pct"].get().strip()
        args.hotspot_pct = float(pct) if pct else None
        args.keep_terminal_gaps = not state["trim"].get()
        args.title = state["title"].get().strip() or None
        args.lang = state["lang"].get()
        args.dpi = 100
        args.figsize = None
        args.top = 10
        args.pdf_with_table = True
        args.aliases = read_alias_text()
        for key, value in collect_display_opts().items():
            setattr(args, key, value)
        return args

    def draw():
        path = state["input"].get().strip()
        if not path or not os.path.exists(path):
            messagebox.showwarning(text("msg_missing_title"), text("msg_missing_body"))
            return
        try:
            args = current_args()
        except ValueError as exc:
            messagebox.showerror(text("msg_param_title"), str(exc))
            return
        draw_btn.configure(state="disabled")
        export_btn.configure(state="disabled")
        status.set(text("status_running"))

        def worker():
            try:
                aln = read_alignment(args.input)
                clashes = duplicate_display_names(aln.names, args.aliases)
                if clashes:
                    task_queue.put(("dup", clashes, None))
                    return
                result = analyse(
                    aln,
                    coords_seq=args.coords_seq,
                    window=args.window,
                    step=args.step,
                    indel_min=args.indel_min,
                    hotspot_k=args.hotspot_k,
                    hotspot_pct=args.hotspot_pct,
                    keep_terminal_gaps=args.keep_terminal_gaps,
                )
                if args.no_gene_track:
                    # Track hidden: skip the annotation read entirely.
                    result.features = []
                    result.feature_source = "none"
                else:
                    attach_features(
                        result,
                        args.annotation,
                        None,
                        None,
                        90,
                        annotation_options={
                            "timeout": args.annotation_timeout,
                            "use_cache": not args.no_annotation_cache,
                            "refresh": args.refresh_annotation,
                        },
                    )
                task_queue.put(("ok", result, args))
            except Exception as exc:  # noqa: BLE001
                err = "".join(traceback.format_exception_only(type(exc), exc)).strip()
                task_queue.put(("error", err, None))

        threading.Thread(target=worker, daemon=True).start()

    def fail_duplicates(names) -> None:
        draw_btn.configure(state="normal")
        if app["result"] is not None:
            export_btn.configure(state="normal")
        message = text("msg_alias_dup_body").format(names=", ".join(names))
        status.set(text("status_failed").format(msg=message))
        messagebox.showerror(text("msg_alias_dup_title"), message)

    def fail(message: str):
        draw_btn.configure(state="normal")
        if app["result"] is not None:
            export_btn.configure(state="normal")
        status.set(text("status_failed").format(msg=message))
        messagebox.showerror(text("msg_error_title"), message)

    def finish(result, args):
        try:
            configure_matplotlib(args.lang)
            fig = build_figure(result, args)
        except Exception as exc:  # noqa: BLE001
            fail(str(exc))
            return
        if app["canvas"] is not None:
            app["canvas"].get_tk_widget().destroy()
            app["canvas"] = None
        placeholder.pack_forget()
        canvas = FigureCanvasTkAgg(fig, master=tab_figure)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        app["canvas"] = canvas
        app["fig"] = fig
        app["result"] = result
        app["args"] = args
        draw_btn.configure(state="normal")
        export_btn.configure(state="normal")
        refresh_table_headers()
        refresh_table_rows()
        refresh_table_layout()
        status.set(
            text("status_done").format(
                nseq=result.alignment.nseq,
                cols=result.core_columns,
                ident=result.identity_ungapped,
                hotspots=len(result.hotspots),
            )
        )

    def poll_tasks():
        try:
            while True:
                kind, payload, extra = task_queue.get_nowait()
                if kind == "ok":
                    finish(payload, extra)
                elif kind == "dup":
                    fail_duplicates(payload)
                else:
                    fail(payload)
        except queue.Empty:
            pass
        finally:
            root.after(120, poll_tasks)

    def do_export():
        result = app["result"]
        args = app["args"]
        fig = app["fig"]
        if result is None or args is None or fig is None:
            messagebox.showwarning(text("msg_export_title"), text("msg_export_no_result"))
            return
        want_pdf = state["exp_pdf_fig"].get() or state["exp_pdf_table"].get()
        selected = any(
            (
                state["exp_data"].get(),
                state["exp_svg"].get(),
                state["exp_png"].get(),
                state["exp_md"].get(),
                want_pdf,
            )
        )
        if not selected:
            messagebox.showwarning(text("msg_export_title"), text("msg_export_none"))
            return
        dpi = None
        if state["exp_png"].get():
            try:
                dpi = int(state["png_dpi"].get().strip())
            except ValueError:
                messagebox.showerror(text("msg_export_title"), text("msg_dpi_invalid"))
                return
            if not 72 <= dpi <= 1200:
                messagebox.showerror(text("msg_export_title"), text("msg_dpi_invalid"))
                return
        folder = filedialog.askdirectory(title=text("msg_dir_title"))
        if not folder:
            return
        stem = os.path.splitext(os.path.basename(state["input"].get()))[0] or "alignment"
        written: List[str] = []
        # Exports honour the display switches as they are now, not as they were
        # when the preview was drawn.
        try:
            fresh = current_args()
        except ValueError as exc:
            messagebox.showerror(text("msg_param_title"), str(exc))
            return
        for key, value in collect_display_opts().items():
            setattr(args, key, value)
        args.aliases = fresh.aliases
        args.title = fresh.title
        # Export from a freshly built figure so the file size does not depend on
        # how the preview canvas happened to be resized.
        try:
            export_fig = build_figure(result, args)
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror(text("msg_export_failed_title"), str(exc))
            return
        try:
            if state["exp_data"].get():
                fmt = state["exp_data_format"].get().upper()
                if fmt == "JSON":
                    data_path = os.path.join(folder, f"{stem}_summary.json")
                    with open(data_path, "w", encoding="utf-8") as handle:
                        json.dump(summary_dict(result, args), handle, ensure_ascii=False, indent=2)
                else:
                    data_path = os.path.join(folder, f"{stem}_windows.tsv")
                    write_windows_tsv(result, data_path)
                written.append(data_path)
            if state["exp_svg"].get():
                svg_path = os.path.join(folder, f"{stem}_conservation.svg")
                export_fig.savefig(svg_path, format="svg", bbox_inches="tight")
                written.append(svg_path)
            if state["exp_png"].get():
                png_path = os.path.join(folder, f"{stem}_conservation.png")
                export_fig.savefig(png_path, format="png", dpi=dpi, bbox_inches="tight")
                written.append(png_path)
            if state["exp_md"].get():
                md_path = os.path.join(folder, f"{stem}_hotspots.md")
                write_hotspot_markdown(result, md_path, args, table_lang())
                written.append(md_path)
            if want_pdf:
                pdf_path = os.path.join(folder, f"{stem}_conservation.pdf")
                save_pdf_bundle(
                    export_fig,
                    result,
                    args,
                    pdf_path,
                    include_figure=state["exp_pdf_fig"].get(),
                    include_table=state["exp_pdf_table"].get(),
                    lang=table_lang(),
                )
                written.append(pdf_path)
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror(text("msg_export_failed_title"), str(exc))
            status.set(text("status_failed").format(msg=str(exc)))
            return
        status.set(text("status_exported").format(count=len(written), path=folder))
        # Every file lands in the same folder, so list bare names and show the
        # folder once instead of repeating the full path for each file.
        names = [os.path.basename(path) for path in written]
        messagebox.showinfo(
            text("msg_export_done_title"),
            text("msg_export_done_body").format(
                files=text("export_list_sep").join(names),
                path=folder,
            ),
        )

    draw_btn.configure(command=draw)
    export_btn.configure(command=do_export)

    def on_close():
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    def measure_panel_width() -> int:
        """Widest needed width across both languages, so switching never reflows."""
        saved = gui_lang["value"]
        widths = []
        for lang_value in ("en", "zh"):
            gui_lang["value"] = lang_value
            for widget, key, attr in registry:
                try:
                    widget.configure(**{attr: text(key)})
                except tk.TclError:
                    pass
            root.update_idletasks()
            widths.append(left.winfo_reqwidth())
        gui_lang["value"] = saved
        return min(max(widths) + 24, 640)

    left_canvas.configure(width=measure_panel_width())

    apply_language()
    sync_gene_dependent_widgets()
    sync_note_and_indel()
    refresh_table_layout()

    def _bind_wheel_tree(widget) -> None:
        """Widget-level bindings fire before class bindings, which may break."""
        for child in widget.winfo_children():
            if child is not alias_text:
                for sequence in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                    child.bind(sequence, _on_wheel, add="+")
                if touchpad_events:
                    child.bind("<TouchpadScroll>", _on_touchpad, add="+")
                _bind_wheel_tree(child)

    _bind_wheel_tree(left)
    left_canvas.yview_moveto(0.0)
    root.after(120, poll_tasks)

    if self_test is not None:
        self_test(
            {
                "root": root,
                "app": app,
                "state": state,
                "tree": tree,
                "notebook": notebook,
                "status": status,
                "current_args": current_args,
                "draw": draw,
                "finish": finish,
                "do_export": do_export,
                "apply_language": apply_language,
                "toggle_language": toggle_language,
                "gui_lang": gui_lang,
                "alias_text": alias_text,
                "fill_aliases": fill_aliases,
                "anno_entry": anno_entry,
                "browse_ann": browse_ann,
                "column_checks": column_checks,
                "display_checks": display_checks,
                "refresh_table_layout": refresh_table_layout,
                "table_title_label": table_title_label,
                "table_subtitle_label": table_subtitle_label,
                "table_note": table_note,
                "left_canvas": left_canvas,
            }
        )
        root.destroy()
        return 0

    root.mainloop()
    return 0


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.show_font:
        configure_matplotlib(args.lang)
        chosen = pick_figure_font(available_font_families())
        print(font_notice(chosen, args.lang))
        return 0
    if args.gui or (args.input is None and args.out is None):
        return launch_gui()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())

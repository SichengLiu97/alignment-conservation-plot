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
        --coords-seq 2 --annotation ref.gb

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

    match_pairs = 0.0
    pair_total = 0.0
    gaps = 0
    for i in core_slice:
        same, total = column_pair_score(seqs, i)
        if total == 0:
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
            if total == 0:
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
        for group in merged:
            if len(group) < 2 and len(low) > 1:
                # keep isolated single windows only when they are clearly worst
                if group[0].ident_ungapped > threshold - 1.0:
                    continue
            m = sum(g.matches for g in group)
            mm = sum(g.mismatches for g in group)
            gp = sum(g.gaps for g in group)
            n = sum(g.ncols for g in group)
            hotspots.append(
                Hotspot(
                    start=group[0].start,
                    end=group[-1].end,
                    ident=100.0 * m / (m + mm) if m + mm else 0.0,
                    gap_frac=100.0 * gp / n if n else 0.0,
                    ncols=n,
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
    from Bio import SeqIO

    handle_ctx = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")
    with handle_ctx as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records:
        return []
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
) -> None:
    aln = result.alignment
    if path:
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
            result.warnings.append(f"注释文件 {os.path.basename(path)} 中没有可用特征")
        return

    seq = aln.ungapped(result.coords_seq)
    if len(seq) >= 300 and set(seq.upper()) <= set("ACGTN"):
        result.features = predict_orfs(seq, min_orf_aa)
        result.feature_source = "predicted ORFs"
    else:
        result.feature_source = "none"


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
        "footer": (
            "窗口 {window} bp / 步长 {step} bp；一致度不计 gap 列，中位数 {median:.1f}%。"
            "已剔除两端全 gap 列 {lead:,} / {trail:,} 列。基因来源：{source}。"
            "比较区间：{start:,}–{end:,}（{span:,} bp）。"
        ),
    },
}


def configure_matplotlib(lang: str) -> None:
    import logging

    import matplotlib

    # PingFang SC has no separate bold face; matplotlib falls back silently but
    # logs a noisy warning on every run.
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

    matplotlib.rcParams["font.sans-serif"] = [
        "PingFang SC",
        "Hiragino Sans GB",
        "Heiti SC",
        "Songti SC",
        "Source Han Sans SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
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
    ref_name = aln.names[result.coords_seq]
    xmin, xmax = result.coord_min, result.coord_max

    n_genes = len(result.features)
    row_heights = [0.62, 3.4, 0.85] + ([0.9] if n_genes else [])
    if opts.figsize:
        figsize = tuple(opts.figsize)
    else:
        figsize = (11.0, 5.7 if n_genes else 4.6)
    fig = Figure(figsize=figsize, dpi=opts.dpi, layout="constrained")

    gs = fig.add_gridspec(
        nrows=len(row_heights),
        ncols=1,
        height_ratios=row_heights,
        hspace=0.10,
    )
    ax_head = fig.add_subplot(gs[0])
    ax_head.axis("off")
    ax = fig.add_subplot(gs[1], sharex=ax_head)
    ax_gap = fig.add_subplot(gs[2], sharex=ax)

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

    # hotspot numbers
    for i, hs in enumerate(result.hotspots[:9], start=1):
        ax.text(
            (hs.start + hs.end) / 2,
            yhi - (yhi - ylo) * 0.06,
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

    # indel ticks
    if result.indels:
        ymax = ax_gap.get_ylim()[1]
        for ind in result.indels:
            ax_gap.plot([ind.coord, ind.coord], [0, ymax * 0.18], color="#6a3d9a", lw=1.0, alpha=0.8)

    # --- gene track --------------------------------------------------------
    ax_gene = None
    if n_genes:
        ax_gene = fig.add_subplot(gs[3], sharex=ax)
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
    bottom_ax = ax_gene if ax_gene is not None else ax_gap
    bottom_ax.set_xlabel(t["pos"].format(name=ref_name), fontsize=9.5)

    def fmt_coord(v, _pos):
        return f"{int(v):,}"

    bottom_ax.xaxis.set_major_formatter(FuncFormatter(fmt_coord))

    # --- titles and footer -------------------------------------------------
    if opts.title:
        title = opts.title
    elif aln.nseq == 2:
        title = f"{aln.names[0]} vs {aln.names[1]}"
    else:
        title = f"{aln.names[0]} (+{aln.nseq - 1} more)"
    subtitle_bits = []
    if aln.nseq > 2:
        subtitle_bits.append(f"{aln.nseq} sequences, mean pairwise identity")
    subtitle_bits.append(f"{result.core_columns:,} alignment columns")
    subtitle_bits.append(f"{result.identity_ungapped:.1f}% identity")
    ax_head.text(0.0, 0.92, title, fontsize=12.5, fontweight="bold", ha="left", va="top")
    ax_head.text(0.0, 0.18, "  |  ".join(subtitle_bits), fontsize=8, color="#555", ha="left", va="bottom")

    footer = t["footer"].format(
        window=opts.window,
        step=opts.step,
        median=result.identity_ungapped,
        lead=result.lead_gap,
        trail=result.trail_gap,
        source=result.feature_source,
        start=xmin,
        end=xmax,
        span=xmax - xmin + 1,
    )
    fig.text(0.01, -0.02, footer, fontsize=7, color="#666", ha="left", va="top")
    return fig


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


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


def summary_dict(result: AnalysisResult, opts) -> Dict[str, object]:
    aln = result.alignment
    return {
        "sequences": aln.names,
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
                "start": h.start,
                "end": h.end,
                "identity_pct": round(h.ident, 2),
                "gap_pct": round(h.gap_frac, 2),
                "columns": h.ncols,
            }
            for h in result.hotspots
        ],
        "annotation_source": result.feature_source,
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
    print("=" * 72)
    print(f"sequences        : {', '.join(aln.names)}")
    print(f"format           : {aln.source_format}")
    print(f"alignment columns: {aln.ncol:,}")
    print(
        f"terminal gaps    : {result.lead_gap:,} (5') / {result.trail_gap:,} (3')"
    )
    print(
        f"compared         : {aln.names[result.coords_seq]} "
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
                f"identity {h.ident:5.1f}%  gap {h.gap_frac:5.1f}%"
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
    p.add_argument("--annotation", help="GenBank / GFF3 / TSV annotation file")
    p.add_argument("--annotation-format", choices=["genbank", "gff3", "tsv"], help="override format detection")
    p.add_argument("--record", help="GenBank record id to use for annotation")
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
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
        print(f"PDF written: {pdf_path}")
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
    aln = read_alignment(args.input)
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
    attach_features(
        result,
        args.annotation,
        args.annotation_format,
        args.record,
        args.min_orf_aa,
    )
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
    return 0


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------


def launch_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    from tkinter import font as tkfont

    configure_matplotlib("zh")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

    root = tk.Tk()
    root.title("Alignment Conservation Plot")
    root.geometry("1280x820")
    try:
        tkfont.nametofont("TkDefaultFont").configure(family="PingFang SC", size=12)
        tkfont.nametofont("TkTextFont").configure(family="PingFang SC", size=12)
        tkfont.nametofont("TkMenuFont").configure(family="PingFang SC", size=12)
    except Exception:
        pass

    state = {
        "input": tk.StringVar(),
        "annotation": tk.StringVar(),
        "coords_seq": tk.StringVar(value="2"),
        "window": tk.StringVar(value="300"),
        "step": tk.StringVar(value="50"),
        "indel_min": tk.StringVar(value="8"),
        "hotspot_k": tk.StringVar(value="1.5"),
        "hotspot_pct": tk.StringVar(value=""),
        "title": tk.StringVar(),
        "lang": tk.StringVar(value="zh"),
        "trim": tk.BooleanVar(value=True),
    }

    left = ttk.Frame(root, padding=(12, 12, 8, 12))
    left.pack(side="left", fill="y")
    right = ttk.Frame(root, padding=(4, 8, 8, 8))
    right.pack(side="right", fill="both", expand=True)

    def labelled_entry(parent, row, text, var, width=28):
        ttk.Label(parent, text=text).grid(row=row, column=0, sticky="w", pady=3)
        entry = ttk.Entry(parent, textvariable=var, width=width)
        entry.grid(row=row, column=1, columnspan=2, sticky="we", pady=3)
        return entry

    row = 0
    ttk.Label(left, text="比对文件", font=("PingFang SC", 12, "bold")).grid(row=row, column=0, columnspan=3, sticky="w")
    row += 1
    entry_in = ttk.Entry(left, textvariable=state["input"], width=34)
    entry_in.grid(row=row, column=0, columnspan=2, sticky="we", pady=3)

    def pick_input():
        path = filedialog.askopenfilename(
            title="选择比对文件",
            filetypes=[
                ("Alignment", "*.srspair *.pair *.aln *.fasta *.fa *.fas *.sto *.phy *.msf *.nex *.gz"),
                ("All files", "*.*"),
            ],
        )
        if path:
            state["input"].set(path)

    ttk.Button(left, text="浏览…", command=pick_input, width=8).grid(row=row, column=2, padx=(4, 0))
    row += 1

    ttk.Label(left, text="注释文件（可空）").grid(row=row, column=0, columnspan=3, sticky="w", pady=(10, 0))
    row += 1
    ttk.Entry(left, textvariable=state["annotation"], width=34).grid(row=row, column=0, columnspan=2, sticky="we", pady=3)

    def pick_annotation():
        path = filedialog.askopenfilename(
            title="选择注释文件",
            filetypes=[("Annotation", "*.gb *.gbk *.genbank *.gff *.gff3 *.tsv *.txt"), ("All files", "*.*")],
        )
        if path:
            state["annotation"].set(path)

    ttk.Button(left, text="浏览…", command=pick_annotation, width=8).grid(row=row, column=2, padx=(4, 0))
    row += 1

    ttk.Separator(left, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="we", pady=10)
    row += 1
    ttk.Label(left, text="参数", font=("PingFang SC", 12, "bold")).grid(row=row, column=0, columnspan=3, sticky="w")
    row += 1

    labelled_entry(left, row, "坐标参考序列 (1-based)", state["coords_seq"], width=12); row += 1
    labelled_entry(left, row, "窗口大小 (bp)", state["window"], width=12); row += 1
    labelled_entry(left, row, "步长 (bp)", state["step"], width=12); row += 1
    labelled_entry(left, row, "indel 阈值 (bp)", state["indel_min"], width=12); row += 1
    labelled_entry(left, row, "热点 k (MAD)", state["hotspot_k"], width=12); row += 1
    labelled_entry(left, row, "热点百分位 (空=用 k)", state["hotspot_pct"], width=12); row += 1
    labelled_entry(left, row, "标题（可空）", state["title"], width=28); row += 1

    ttk.Checkbutton(left, text="剔除两端全 gap 列", variable=state["trim"]).grid(
        row=row, column=0, columnspan=3, sticky="w", pady=3
    )
    row += 1
    lang_frame = ttk.Frame(left)
    lang_frame.grid(row=row, column=0, columnspan=3, sticky="w", pady=3)
    ttk.Label(lang_frame, text="图内文字：").pack(side="left")
    ttk.Radiobutton(lang_frame, text="中文", value="zh", variable=state["lang"]).pack(side="left")
    ttk.Radiobutton(lang_frame, text="English", value="en", variable=state["lang"]).pack(side="left")
    row += 1

    ttk.Separator(left, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="we", pady=10)
    row += 1
    status = tk.StringVar(value="就绪")
    buttons = ttk.Frame(left)
    buttons.grid(row=row, column=0, columnspan=3, sticky="we", pady=(2, 6))
    draw_btn = ttk.Button(buttons, text="绘图", width=10)
    draw_btn.pack(side="left")
    export_svg = ttk.Button(buttons, text="导出 SVG", width=11, state="disabled")
    export_svg.pack(side="left", padx=4)
    export_pdf = ttk.Button(buttons, text="导出 PDF", width=11, state="disabled")
    export_pdf.pack(side="left", padx=4)
    export_data = ttk.Button(buttons, text="导出数据", width=11, state="disabled")
    export_data.pack(side="left", padx=4)
    row += 1
    ttk.Label(left, textvariable=status, foreground="#555", wraplength=340).grid(
        row=row, column=0, columnspan=3, sticky="w"
    )

    holder = ttk.Frame(right)
    holder.pack(fill="both", expand=True)
    placeholder = ttk.Label(holder, text="请选择比对文件后点击「绘图」", foreground="#888")
    placeholder.pack(expand=True)

    app = {"fig": None, "result": None, "canvas": None, "toolbar": None}

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
        return args

    def draw():
        path = state["input"].get().strip()
        if not path or not os.path.exists(path):
            messagebox.showwarning("缺少输入", "请先选择一个存在的比对文件。")
            return
        try:
            args = current_args()
        except ValueError as exc:
            messagebox.showerror("参数错误", str(exc))
            return
        draw_btn.configure(state="disabled")
        status.set("正在解析并计算…")

        def worker():
            try:
                aln = read_alignment(args.input)
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
                attach_features(result, args.annotation, None, None, 90)
                root.after(0, lambda: finish(result, args, None))
            except Exception as exc:  # noqa: BLE001
                err = "".join(traceback.format_exception_only(type(exc), exc)).strip()
                root.after(0, lambda: fail(err))

        import threading

        threading.Thread(target=worker, daemon=True).start()

    def fail(message: str):
        draw_btn.configure(state="normal")
        status.set("失败：" + message)
        messagebox.showerror("绘图失败", message)

    def finish(result, args, _unused):
        try:
            configure_matplotlib(args.lang)
            fig = build_figure(result, args)
        except Exception as exc:  # noqa: BLE001
            fail(str(exc))
            return
        if app["canvas"] is not None:
            app["canvas"].get_tk_widget().destroy()
        if app["toolbar"] is not None:
            app["toolbar"].destroy()
            app["toolbar"] = None
        placeholder.pack_forget()
        canvas = FigureCanvasTkAgg(fig, master=holder)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        app["canvas"] = canvas
        app["fig"] = fig
        app["result"] = result
        draw_btn.configure(state="normal")
        for btn in (export_svg, export_pdf, export_data):
            btn.configure(state="normal")
        status.set(
            f"完成：{result.alignment.nseq} 条序列，{result.core_columns:,} 列，"
            f"一致度 {result.identity_ungapped:.1f}%，热点 {len(result.hotspots)} 处"
        )

    def save_fig(fmt: str):
        if app["fig"] is None:
            return
        default = os.path.splitext(os.path.basename(state["input"].get()))[0] + "_conservation." + fmt
        path = filedialog.asksaveasfilename(
            title=f"导出 {fmt.upper()}",
            defaultextension="." + fmt,
            initialfile=default,
            filetypes=[(fmt.upper(), "*." + fmt)],
        )
        if not path:
            return
        try:
            app["fig"].savefig(path, format=fmt, bbox_inches="tight", dpi=300)
            status.set(f"已导出：{path}")
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("导出失败", str(exc))

    def save_data():
        if app["result"] is None:
            return
        base = os.path.splitext(os.path.basename(state["input"].get()))[0]
        path = filedialog.asksaveasfilename(
            title="导出窗口统计",
            defaultextension=".tsv",
            initialfile=base + "_windows.tsv",
            filetypes=[("TSV", "*.tsv"), ("JSON", "*.json")],
        )
        if not path:
            return
        try:
            args = current_args()
            if path.lower().endswith(".json"):
                with open(path, "w", encoding="utf-8") as handle:
                    json.dump(summary_dict(app["result"], args), handle, ensure_ascii=False, indent=2)
            else:
                write_windows_tsv(app["result"], path)
            status.set(f"已导出：{path}")
        except Exception as exc:  # noqa: BLE001
            messagebox.showerror("导出失败", str(exc))

    draw_btn.configure(command=draw)
    export_svg.configure(command=lambda: save_fig("svg"))
    export_pdf.configure(command=lambda: save_fig("pdf"))
    export_data.configure(command=save_data)

    def on_close():
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.mainloop()
    return 0


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.gui or (args.input is None and args.out is None):
        return launch_gui()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())

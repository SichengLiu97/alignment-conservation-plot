# Example output

[中文](README_zh.md) | **English**

## Provenance, acknowledgement and licence

**Data source and toolchain.** Every input is public, and every step ran on
[Galaxy Australia](https://usegalaxy.org.au):

- **Reads** — NCBI SRA run [`SRR11528307`](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR11528307&display=metadata) (BioProject `PRJNA624792`):
  genome sequencing of a SARS-CoV-2 clinical isolate (SB2, Canada) after one
  passage in Vero E6 cells.
- **Assembly** — `Assembly / SPAdes` (Galaxy), producing the 20,761 bp contig
  `NODE_1_length_20761_cov_2834.337024`. SPAdes is not fully deterministic, so
  a re-assembly may not reproduce this contig exactly.
- **Reference** — `Get Data / NCBI Accession Download`: GenBank `NC_004718.3`
  (SARS coronavirus Tor2, complete genome).
- **Alignment** — `EMBOSS / needle` with `-gapopen 10.0 -gapextend 0.5
  -aformat3 srspair`.

**Acknowledgement.** When you reuse this example, please acknowledge Galaxy
Australia:

> "This work is supported by Galaxy Australia, a service provided by Australian
> BioCommons and its partners. The service receives NCRIS funding through
> Bioplatforms Australia, as well as The University of Melbourne and
> Queensland Government RICF funding."

and cite the platform paper:

> The Galaxy Community. *Galaxy for accessible, reproducible, and collaborative
> data analyses: 2026 update*. Nucleic Acids Research, 54(W1):W105–W116, 2026.
> <https://doi.org/10.1093/nar/gkag469>

**Licence.** The BSD-3-Clause licence covering this repository applies to the
**source code only**. Nothing in this directory — neither the `.srspair`
alignment nor the generated `example.svg`, `example.pdf`, `example.png`,
`windows.tsv`, `summary.json` and `hotspots.md` — is licensed under it. Those
files are derived from the public data listed above: NCBI places no
restrictions on the use or distribution of the data it serves, and INSDC
records carry no use restrictions or licensing requirements. Cite the original
data sources if you redistribute them.

## Files and how to regenerate them

Everything in this directory was generated from
`examples/sars-cov-2-assembled-contig-vs-sars-cov.srspair`:

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "examples/sars-cov-2-assembled-contig-vs-sars-cov.srspair" \
    -o examples/example.svg --pdf --png \
    --coords-seq 2 --annotation NC_004718.3 \
    --alias "NODE_1_length=SRR11528307 Contig 1" \
    --stats-out examples/windows.tsv --summary-out examples/summary.json \
    --hotspot-table-md examples/hotspots.md --pdf-with-table
```

| File | Description |
|---|---|
| `example.svg` / `example.pdf` / `example.png` | The same figure as SVG / PDF / PNG |
| `windows.tsv` | Per-window statistics |
| `summary.json` | Summary (sequences, terminal gaps, hotspots, indels, genes) |
| `hotspots.md` | Low-conservation table as Markdown (rebuild with `--hotspot-table-md`) |
| `sars-cov-2-assembled-contig-vs-sars-cov.srspair` | Example input: a needle alignment between the SARS-CoV reference genome and one SARS-CoV-2 contig |

The annotation file is not part of the repository; pass the NCBI accession
instead and the program downloads it:

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "examples/sars-cov-2-assembled-contig-vs-sars-cov.srspair" \
    -o examples/example.svg --pdf --coords-seq 2 \
    --annotation NC_004718.3 \
    --alias "NODE_1_length=SRR11528307 Contig 1"
```

The first run caches the record at
`~/Library/Caches/align-conservation-plot/NC_004718.3.gb`, so later runs work
offline.

The `--alias` flag renames the contig in the figure only; the alignment file is
untouched, and the data exports keep the original name alongside the alias.

The figures are drawn with the first open-licence font available on your system
(here "Source Han Sans SC"); `--show-font` reports which one and its licence —
see [Fonts](../README.md#fonts).

Numbers to compare against: terminal gaps 8,917 / 119, 21,072 core columns,
17,219 matches, 3,185 mismatches, identity 84.4%.

> Binary files (pdf/png) are excluded from version control by repository
> convention; rebuild them with the commands above whenever needed.

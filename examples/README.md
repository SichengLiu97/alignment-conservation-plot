# Example output

[中文](README_zh.md) | **English**

Generated from `examples/sars-cov-vs-sars-cov-2.srspair`:

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "examples/sars-cov-vs-sars-cov-2.srspair" \
    -o examples/example.svg --pdf --png \
    --coords-seq 2 --annotation NC_004718.3 \
    --stats-out examples/windows.tsv --summary-out examples/summary.json \
    --hotspot-table-md examples/hotspots.md --pdf-with-table
```

| File | Description |
|---|---|
| `example.svg` / `example.pdf` / `example.png` | The same figure as SVG / PDF / PNG |
| `windows.tsv` | Per-window statistics |
| `summary.json` | Summary (sequences, terminal gaps, hotspots, indels, genes) |
| `hotspots.md` | Low-conservation table as Markdown (rebuild with `--hotspot-table-md`) |
| `sars-cov-vs-sars-cov-2.srspair` | Example input: a needle alignment between the SARS-CoV reference genome and one SARS-CoV-2 contig |

The annotation file is not part of the repository; pass the NCBI accession
instead and the program downloads it:

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "examples/sars-cov-vs-sars-cov-2.srspair" \
    -o examples/example.svg --pdf --coords-seq 2 \
    --annotation NC_004718.3
```

The first run caches the record at
`~/Library/Caches/align-conservation-plot/NC_004718.3.gb`, so later runs work
offline.

Numbers to compare against: terminal gaps 8,917 / 119, 21,072 core columns,
17,219 matches, 3,185 mismatches, identity 84.4%.

> Binary files (pdf/png) are excluded from version control by repository
> convention; rebuild them with the commands above whenever needed.

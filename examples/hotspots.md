# Low-conservation regions

- Alignment: `NODE_1_length` vs `NC_004718.3`
- Coordinate sequence: `NC_004718.3`
- Compared region: 8,918–29,632 (20,715 bp)
- Window / step: 300 bp / 50 bp
- Overall identity: 84.4% (gap-free)
- Gene source: GenBank

| # | Coordinates (bp) | Length (bp) | Identity % (ungapped) | Identity % (with gaps) | Indel % (gap cols / indel bp) | Region |
|---|---|---|---|---|---|---|
| 1 | 9,168–9,617 | 450 | 77.10 | 74.24 | 3.7 / 0.0 | ORF1ab |
| 2 | 21,518–23,817 | 2300 | 75.69 | 67.57 | 10.7 / 5.4 | S |
| 3 | 25,318–26,017 | 700 | 76.71 | 73.60 | 4.1 / 0.0 | ORF3a / ORF3b |

Indel % shows gap-column share / share covered by indels at or above the reporting threshold.

```bash
python align_conservation_plot.py -i <alignment> -o <figure.svg> --coords-seq 2 --window 300 --step 50
```

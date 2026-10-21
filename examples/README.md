# 示例输出

用 `../../1 - Assembly and Annotation of the SARS-CoV-2 Genome/Galaxy17-[needle on dataset 9 and 16].srspair`
生成：

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "../../1 - Assembly and Annotation of the SARS-CoV-2 Genome/Galaxy17-[needle on dataset 9 and 16].srspair" \
    -o examples/example.svg --pdf --png \
    --coords-seq 2 --annotation examples/reference.gb \
    --stats-out examples/windows.tsv --summary-out examples/summary.json \
    --hotspot-table-md examples/hotspots.md --pdf-with-table
```

| 文件 | 说明 |
|---|---|
| `example.svg` / `example.pdf` / `example.png` | 同一张图，SVG / PDF / PNG |
| `windows.tsv` | 逐窗口统计 |
| `summary.json` | 汇总（序列、端部 gap、热点、插入或删除、基因） |
| `hotspots.md` | 低保守区表（Markdown，可用 `--hotspot-table-md` 重建） |
| `reference.gb` | NC_004718.3 的 GenBank 注释，供基因轨道使用（离线复现用） |

也可以不下载注释文件，直接给 accession：

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "../../1 - Assembly and Annotation of the SARS-CoV-2 Genome/Galaxy17-[needle on dataset 9 and 16].srspair" \
    -o examples/example.svg --pdf --coords-seq 2 \
    --annotation NC_004718.3
```

首次运行会从 NCBI 下载并缓存到 `~/Library/Caches/align-conservation-plot/NC_004718.3.gb`，
之后离线也能用。

核对数值：端部 gap 8,917 / 119，核心区 21,072 列，matches 17,219，
mismatches 3,185，identity 84.4%。

> 二进制文件（pdf/png）按仓库约定不纳入版本控制，可用上面的命令随时重建。

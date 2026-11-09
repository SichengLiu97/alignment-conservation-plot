# 示例输出

**中文** | [English](README.md)

## 数据来源、致谢与许可

**数据来源与工具链。** 全部输入均为公开数据，全部步骤都在
[Galaxy Australia](https://usegalaxy.org.au) 上完成：

- **测序数据** —— NCBI SRA run [`SRR11528307`](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR11528307&display=metadata)
  （BioProject `PRJNA624792`）：一株 SARS-CoV-2 临床分离株（SB2，加拿大）
  在 Vero E6 细胞中传代一次后的基因组测序数据。
- **组装** —— `Assembly / SPAdes`（Galaxy），得到 20,761 bp 的 contig
  `NODE_1_length_20761_cov_2834.337024`。SPAdes 并非完全确定性算法，
  重新组装未必能复现同一条 contig。
- **参考序列** —— `Get Data / NCBI Accession Download`：GenBank
  `NC_004718.3`（SARS coronavirus Tor2，完整基因组）。
- **比对** —— `EMBOSS / needle`，参数 `-gapopen 10.0 -gapextend 0.5
  -aformat3 srspair`。

**致谢。** 复用本示例时，请致谢 Galaxy Australia：

> 「This work is supported by Galaxy Australia, a service provided by Australian
> BioCommons and its partners. The service receives NCRIS funding through
> Bioplatforms Australia, as well as The University of Melbourne and
> Queensland Government RICF funding.」

并引用其平台论文：

> The Galaxy Community. *Galaxy for accessible, reproducible, and collaborative
> data analyses: 2026 update*. Nucleic Acids Research, 54(W1):W105–W116, 2026.
> <https://doi.org/10.1093/nar/gkag469>

**许可。** 本仓库的 BSD-3-Clause 许可证仅适用于**源代码**。本目录下的文件
——包括 `.srspair` 比对文件，以及生成的 `example.svg`、`example.pdf`、
`example.png`、`windows.tsv`、`summary.json`、`hotspots.md`——都**不**在该
许可证的覆盖范围内。它们派生自上述公开数据：NCBI 对其提供的数据的使用与再
分发不设限制，INSDC 记录也不附带任何使用限制或许可要求。若要再分发这些
文件，请引用原始数据来源。

用 `examples/sars-cov-2-assembled-contig-vs-sars-cov.srspair`
生成：

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "examples/sars-cov-2-assembled-contig-vs-sars-cov.srspair" \
    -o examples/example.svg --pdf --png \
    --coords-seq 2 --annotation NC_004718.3 \
    --stats-out examples/windows.tsv --summary-out examples/summary.json \
    --hotspot-table-md examples/hotspots.md --pdf-with-table
```

| 文件 | 说明 |
|---|---|
| `example.svg` / `example.pdf` / `example.png` | 同一张图，SVG / PDF / PNG |
| `windows.tsv` | 逐窗口统计 |
| `summary.json` | 汇总（序列、端部 gap、热点、插入或删除、基因） |
| `hotspots.md` | 低保守区表（Markdown，可用 `--hotspot-table-md` 重建） |
| `sars-cov-2-assembled-contig-vs-sars-cov.srspair` | 示例输入：SARS-CoV 参考基因组与一条 SARS-CoV-2 contig 的 needle 比对 |

注释文件不进仓库，直接给 NCBI accession 即可（程序会自动下载）：

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "examples/sars-cov-2-assembled-contig-vs-sars-cov.srspair" \
    -o examples/example.svg --pdf --coords-seq 2 \
    --annotation NC_004718.3
```

首次运行会从 NCBI 下载并缓存到 `~/Library/Caches/align-conservation-plot/NC_004718.3.gb`，
之后离线也能用。

核对数值：端部 gap 8,917 / 119，核心区 21,072 列，matches 17,219，
mismatches 3,185，identity 84.4%。

> 二进制文件（pdf/png）按仓库约定不纳入版本控制，可用上面的命令随时重建。

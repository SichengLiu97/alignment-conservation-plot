# align-conservation-plot

从比对文件生成「滑窗一致度 + gap 密度 + 基因轨道」的出版级图，输出 **SVG / PDF / PNG**。
支持 EMBOSS `srspair` / `pair` 以及常见多序列比对格式，带 Tkinter GUI。

![example](examples/example.svg)

## 为什么自己写

现成工具都差一块：

| 工具 | 能力 | 缺口 |
|---|---|---|
| EMBOSS `plotcon` | 滑窗保守度曲线，可输出 pdf/svg | 画的是替换矩阵相似度分数，不是 identity%；无 gap 轨道、无 indel 标记、无基因轨道 |
| MSAexplorer | 多序列分析，读 gb/gff/bed，导 PDF | 无「按参考坐标的滑窗 identity + gap 密度」曲线；仅 PDF |
| pyMSAviz / MSAplot | MSA 色块图 + seqlogo | 不做滑窗 identity 曲线 |
| Jalview | 逐列保守度着色 | 不是按基因组坐标的滑窗图 |
| dna_features_viewer | 基因轨道 | 仅轨道 |

`srspair` 格式目前没有专门的可视化工具（GitHub 全站检索只在 janis-core 里找到 datatype 声明）。

## 安装

```bash
cd tools/alignment-conservation-plot
./setup.sh
```

脚本会优先使用 `/opt/homebrew/bin/python3.14`（自带 Tk 9.1），创建 `.venv` 并安装
`matplotlib` 与 `biopython`。若需要代理：

```bash
http_proxy=http://127.0.0.1:8080 https_proxy=http://127.0.0.1:8080 ./setup.sh
```

## 命令行

```bash
./.venv/bin/python align_conservation_plot.py \
    -i "../../1 - Assembly and Annotation of the SARS-CoV-2 Genome/Galaxy17-[needle on dataset 9 and 16].srspair" \
    -o figure.svg --pdf \
    --coords-seq 2 --window 300 --step 50 \
    --annotation NC_004718.3 \
    --title "NODE_1 vs NC_004718.3"
```

输出：

```
figure.svg
figure.pdf
```

### 常用参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `-i, --input` | 必填 | 比对文件（`.gz` 亦可） |
| `-o, --out` | `<input>_conservation.svg` | SVG 输出路径 |
| `--pdf [PATH]` | 关 | 同时输出 PDF |
| `--png [PATH]` | 关 | 同时输出 PNG（`--dpi` 控制分辨率） |
| `--coords-seq` | `1` | 用第几条序列做坐标轴（1-based）。与参考基因组比对时通常设为 `2` |
| `--window` / `--step` | `300` / `50` | 滑窗大小与步长（bp） |
| `--annotation` | 无 | 注释文件路径，**或 NCBI accession（如 `NC_004718.3`）** |
| `--annotation-timeout` | `60` | 在线下载超时（秒） |
| `--no-annotation-cache` | 关 | 不读也不写本地下载缓存 |
| `--refresh-annotation` | 关 | 忽略缓存，重新下载 |
| `--annotation-format` | 自动 | `genbank` / `gff3` / `tsv` |
| `--record` | 第一个 | GenBank 多记录时指定 record id |
| `--min-orf-aa` | `90` | 自动预测 ORF 的最短氨基酸长度 |
| `--indel-min` | `8` | 只报告 ≥ 该长度的 indel |
| `--hotspot-k` | `1.5` | 热点阈值 = median − k × 1.4826 × MAD |
| `--hotspot-pct` | 空 | 改用百分位阈值（优先于 `--hotspot-k`） |
| `--keep-terminal-gaps` | 关 | 保留两端全 gap 列（默认剔除） |
| `--lang` | `en` | 图内文字语言，`zh` 自动用 PingFang SC |
| `--stats-out` | 无 | 导出逐窗口统计 TSV |
| `--summary-out` | 无 | 导出 JSON 摘要 |
| `--figsize W H` | 自动 | 图幅尺寸（英寸） |
| `--gui` | — | 启动图形界面 |

## 图形界面

```bash
./.venv/bin/python align_conservation_plot.py --gui
# 或双击 run_gui.command
```

- 左侧填写参数，右侧实时预览；
- 「注释文件或 accession」框可以直接填 `NC_004718.3`，联网下载（有缓存）；
- 「导出 SVG / PDF / 数据」按钮直接保存；
- 计算在后台线程执行，不阻塞界面。

## 统计口径

- **identity（不计 gap）** = matches / (matches + mismatches)，只在双方都有碱基的列上计算。
- **identity（计 gap）** = matches / 所有列，gap 列计入分母。
- **多序列**：逐列统计「相同的序列对 / 有信息的序列对」，等价于平均两两 identity。
- **端部 gap**：所有序列都为 gap 的连续列，默认剔除并在页脚注明长度。
- **indel**：某条序列中连续 gap ≥ `--indel-min` 的区段。
- **热点**：滑窗 identity 低于 `median − k×1.4826×MAD`（或指定百分位）的窗口，合并为区间。

## 在线注释（accession）

`--annotation` 除了本地文件，也可以直接填 NCBI accession，程序会自动下载：

```bash
./.venv/bin/python align_conservation_plot.py \
    -i aln.srspair -o fig.svg --coords-seq 2 \
    --annotation NC_004718.3
```

- 走 E-utilities `efetch`（`db=nuccore`, `rettype=gbwithparts`）。
- 下载结果缓存在 `~/Library/Caches/align-conservation-plot/<accession>.gb`
  （Linux `$XDG_CACHE_HOME`，Windows `%LOCALAPPDATA%`；可用环境变量
  `ALIGN_PLOT_CACHE` 覆盖）。用 `--refresh-annotation` 强制重下。
- 支持 `NCBI_API_KEY` / `NCBI_EMAIL` 环境变量（提高限速配额）。
- 需要代理时，`http_proxy` / `https_proxy` 会被自动沿用：

  ```bash
  http_proxy=http://127.0.0.1:8080 https_proxy=http://127.0.0.1:8080 \
    ./.venv/bin/python align_conservation_plot.py ... --annotation NC_004718.3
  ```

- 图例页脚会标注来源，例如 `GenBank NC_004718.3 (已下载)` 或 `(缓存)`。

## 注释文件

- **GenBank**（`.gb/.gbk/.genbank`）：取 `gene` / `CDS` / `mat_peptide` 等特征。
- **GFF3**（`.gff/.gff3`）：取 `gene` / `CDS` / `mRNA` 等，坐标按 1-based 闭区间解析。
- **TSV**：`name<TAB>start<TAB>end[<TAB>strand][<TAB>color]`，`#` 开头为注释。
- 不提供注释时，程序会对坐标参考序列做六框 ORF 预测（ATG 起始、贪心去重叠），
  并在页脚标注来源为 `predicted ORFs`。

## 输出示例

`examples/` 下有本仓库自带 srspair 的示例输出（见 `examples/README.md`）。

## 已知限制

- `needle` 默认 `-gapextend 0.5` 会在高分歧区过度调用 1–2 bp 的小 gap；
  这些短 indel 不宜当作生物学信号，建议只看 ≥8 bp 的事件。
- 单条 contig 的低保守区也可能来自组装错误，最好用 reads 回贴验证覆盖度。
- 多序列的 identity 语义是「平均两两」，与「相对第一条序列」不同，图注已说明。
- GUI 不做缩放/平移，改参数后重新绘图即可。

## 目录

```
align_conservation_plot.py   # 全部逻辑：解析 / 统计 / 绘图 / CLI / GUI
requirements.txt             # matplotlib, biopython
setup.sh                     # 建 .venv 并安装依赖
run_gui.command              # macOS 双击启动 GUI
examples/                    # 示例输出
```

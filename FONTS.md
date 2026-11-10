# Fonts and their licences

The figure renderer uses whichever font it finds first in a fixed list. This
page records why that order looks the way it does, and what each font's licence
allows — because a chart destined for publication should not silently depend on
a font that cannot legally travel with it.

Run `./.venv/bin/python align_conservation_plot.py --show-font` to see which
font your installation will actually use. The GUI shows the same line above its
status bar.

## Selection order

Three rules, applied in this order:

1. **Open licence first.** A redistributable font keeps the project's outputs
   reproducible on machines that do not own a commercial font.
2. **Sans-serif before serif.** Screens and journal figures favour sans-serif;
   serif is the last resort.
3. **CJK-capable before Latin-only.** Chinese labels are part of the output, so
   a font that can draw them is preferred even when the figure is in English.

The resulting order:

| # | Font | Script coverage | Licence | Redistributable |
|---|---|---|---|---|
| 1 | Source Han Sans SC | CJK + Latin | SIL OFL 1.1 | yes |
| 2 | Source Han Sans | CJK + Latin | SIL OFL 1.1 | yes |
| 3 | Noto Sans CJK SC | CJK + Latin | SIL OFL 1.1 | yes |
| 4 | Noto Sans CJK | CJK + Latin | SIL OFL 1.1 | yes |
| 5 | Noto Sans | Latin | SIL OFL 1.1 | yes |
| 6 | Source Sans Pro / Source Sans 3 | Latin | SIL OFL 1.1 | yes |
| 7 | DejaVu Sans | Latin | Bitstream Vera licence | yes |
| 8 | PingFang SC / HK / TC | CJK + Latin | proprietary (Apple / DynaComware) | no |
| 9 | Hiragino Sans GB | CJK + Latin | proprietary (Apple / SCREEN) | no |
| 10 | Heiti SC | CJK + Latin | proprietary (Apple / SinoType) | no |
| 11 | Microsoft YaHei / SimHei | CJK + Latin | proprietary (Microsoft) | no |
| 12 | Arial Unicode MS | CJK + Latin | proprietary (Monotype) | no |
| 13 | Songti SC | CJK + Latin | proprietary (Apple / SinoType) | no |
| 14 | Source Han Serif SC | CJK + Latin | SIL OFL 1.1 | chosen only if nothing else exists |

The list is applied *per glyph*: if the chosen font lacks a character,
matplotlib falls back to the next family that has it. That is why a Latin-only
font early in the list still renders Chinese correctly.

## Evidence

Licences were read from the metadata embedded in the installed font files
(`name` table IDs 13 and 14), not from memory:

| Font | Embedded statement |
|---|---|
| Source Han Sans SC | “© 2014-2021 Adobe, with Reserved Font Name 'Source'… licensed under the SIL Open Font License, Version 1.1” |
| PingFang SC | “Copyright © 2015 DynaComware. All rights reserved.” (no licence grant) |
| Hiragino Sans GB | “Copyright © 2007-2023 SCREEN Graphic Solutions Co., Ltd. All Rights Reserved.” |
| Arial Unicode MS | Monotype licence agreement, distributed with Microsoft Office |
| DejaVu Sans | Bitstream Vera Fonts Copyright; “DejaVu changes are in public domain” |
| Heiti SC / Songti SC | “Changzhou SinoType Technology Co., Ltd. All rights reserved.” |

Upstream sources:

- Source Han Sans — <https://github.com/adobe-fonts/source-han-sans>, licence
  text at `LICENSE.txt` (SIL OFL 1.1)
- Noto Sans CJK is the same typeface released under the Noto name by Google,
  also under SIL OFL 1.1
- SIL OFL 1.1 — <https://scripts.sil.org/OFL>
- DejaVu licence — <https://dejavu-fonts.github.io/License.html>

## What this means for outputs

- **SVG** keeps text as text (`svg.fonttype = none`) and records the whole font
  list, so the reader's system resolves it. No font data is redistributed.
- **PDF** embeds a subset of the used font (`pdf.fonttype = 42`). OFL fonts
  allow this. For proprietary system fonts, embedding a subset in a document is
  normally permitted by the operating-system licence, but the font *file* must
  never be redistributed — so this repository does not ship any font.
- **PNG** rasterises glyphs; the same caveat applies as for PDF.
- **Data (`.tsv` / `.json`) and the hotspot Markdown** contain no rendered text
  and are unaffected by the font choice.

## Caveat

Font licences change, and some vendors ship slightly different terms per
region or version. The information here reflects what was installed on the
development machine. Verify it yourself before publishing figures — especially
if you rely on a proprietary font (rows 8–13) rather than an OFL one.

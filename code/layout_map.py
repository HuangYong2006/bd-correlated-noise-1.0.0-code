r"""
Deep layout diagnostic for the compiled PDF.

Previous passes only looked at text spans.  Two things were missed:

  * TABLE RULES.  A tabular's horizontal rules are drawn as vector paths, and a
    table can overflow even when its cell text wraps.  We therefore read
    page.get_drawings() and report the extent of every rule.
  * DOCUMENT ORDER.  Where the appendix headings, the bibliography and each
    float actually land, which is what "the appendix is in the wrong place"
    usually means.

Output:
  1. a page-by-page heading map (detected from font size)
  2. every horizontal rule wider than the column or past the margin
  3. where each figure/table lands relative to the text that references it
  4. the order of the top-level headings

Usage: python code/layout_map.py manuscript/overleaf_bd_paper.pdf
"""

from __future__ import annotations

import re
import sys

import fitz

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

path = sys.argv[1] if len(sys.argv) > 1 else "manuscript/overleaf_bd_paper.pdf"
doc = fitz.open(path)
W = doc[0].rect.width
LEFT, MID, RIGHT = 54.0, 308.0, 562.0

print(f"{path}: {len(doc)} pages, page width {W:.0f} pt, "
      f"text area x=[{LEFT:.0f},{RIGHT:.0f}]")
print("=" * 100)

# ---------------------------------------------------------------- 1. headings
print("\n[1] HEADING MAP  (font size >= 9.5 pt, bold-ish, left of MID or spanning)")
print("-" * 100)
headings = []
for pi, pg in enumerate(doc):
    d = pg.get_text("dict")
    for blk in d.get("blocks", []):
        if blk.get("type") != 0:
            continue
        for line in blk.get("lines", []):
            for s in line.get("spans", []):
                t = s["text"].strip()
                if not t or len(t) < 3:
                    continue
                size = s["size"]
                x0, y0, x1, y1 = s["bbox"]
                if size >= 9.4 and (s["flags"] & 2 ** 4 or size >= 10.5):
                    if 40 < x0 < MID + 60 or x0 < MID:
                        headings.append((pi + 1, round(y0, 1), round(size, 1),
                                         round(x0, 1), t[:74]))
seen = set()
for pi, y, sz, x0, t in headings:
    key = (t,)
    if key in seen:
        continue
    seen.add(key)
    print(f"  p{pi:3d}  y={y:6.1f}  size={sz:5.1f}  x0={x0:6.1f}  {t}")

# ------------------------------------------------------------------ 2. rules
print("\n[2] HORIZONTAL RULES (tabular borders) OUTSIDE THE MARGINS")
print("-" * 100)
n_bad = 0
for pi, pg in enumerate(doc):
    for dr in pg.get_drawings():
        r = dr["rect"]
        w, h = r.width, r.height
        if w < 20 or h > 3:          # not a horizontal rule
            continue
        over = max(r.x1 - RIGHT, LEFT - r.x0)
        if over > 2:
            n_bad += 1
            print(f"  p{pi+1:3d}  rule x=[{r.x0:6.1f},{r.x1:6.1f}] "
                  f"y={r.y0:6.1f}  width={w:6.1f}  past margin by {over:5.1f} pt")
if not n_bad:
    print("  none")
print(f"  total offending rules: {n_bad}")

print("\n[2b] WIDEST RULES PER PAGE (to see how wide each table really is)")
print("-" * 100)
for pi, pg in enumerate(doc):
    rules = [dr["rect"] for dr in pg.get_drawings()
             if dr["rect"].width >= 40 and dr["rect"].height <= 3]
    if rules:
        widest = max(rules, key=lambda r: r.width)
        kind = ("column-width" if widest.width < 300 else
                "full-width" if widest.width < 520 else "TOO WIDE")
        print(f"  p{pi+1:3d}: {len(rules):3d} rules, widest {widest.width:6.1f} pt "
              f"x=[{widest.x0:6.1f},{widest.x1:6.1f}]  -> {kind}")

# --------------------------------------------------------- 3. document order
print("\n[3] DOCUMENT ORDER  (headings + float captions in page order)")
print("-" * 100)
marks = []
for pi, pg in enumerate(doc):
    txt = pg.get_text()
    for pat, lab in ((r"^Appendix [A-Z]", "APPENDIX HEADING"),
                     (r"\bAPPENDIX\b", "APPENDIX (small caps)"),
                     (r"^REFERENCES", "REFERENCES"),
                     (r"^ACKNOWLEDGMENTS", "ACKNOWLEDGMENTS"),
                     (r"^Figure \d+\.", "FIGURE CAPTION"),
                     (r"^TABLE \d+\.", "TABLE CAPTION")):
        for m in re.finditer(pat, txt, re.M):
            ln = txt[m.start():txt.find("\n", m.start())][:70]
            marks.append((pi + 1, lab, ln))
for pi, lab, ln in sorted(marks):
    print(f"  p{pi:3d}  {lab:22s} {ln}")

# ------------------------------------------------------ 4. isolated headings
print("\n[4] PAGE OCCUPANCY  (a page that is mostly empty usually means a float")
print("    was pushed, or an appendix started after a forced page break)")
print("-" * 100)
for pi, pg in enumerate(doc):
    blocks = [b for b in pg.get_text("blocks") if b[4].strip()]
    if not blocks:
        print(f"  p{pi+1:3d}: EMPTY")
        continue
    y1 = max(b[3] for b in blocks)
    frac = y1 / pg.rect.height
    flag = "  <-- mostly empty" if frac < 0.55 else ""
    print(f"  p{pi+1:3d}: text ends at y={y1:6.1f} ({frac*100:4.0f}% of page){flag}")

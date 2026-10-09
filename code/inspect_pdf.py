r"""
Diagnose the layout of the compiled PDF.

We cannot display images, so we extract geometry and text instead:
  * page count and page size
  * where the text blocks live, and which ones stick out of the text area
    (a proxy for Overfull \hbox, which is what "the format looks wrong" usually
    means)
  * the number of columns on each page
  * width of every table/figure region
  * unresolved references ("??"), unresolved citations ("[?]") and leftover
    TODO markers
  * the header/footer, and the first lines of the document

Usage: python code/inspect_pdf.py manuscript/overleaf_bd_paper.pdf
"""

from __future__ import annotations

import re
import sys

import fitz

# the Windows console here is GBK; force UTF-8 so Unicode maths does not crash
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

path = sys.argv[1] if len(sys.argv) > 1 else "manuscript/overleaf_bd_paper.pdf"
doc = fitz.open(path)

print(f"file       : {path}")
print(f"pages      : {len(doc)}")
p0 = doc[0]
print(f"page size  : {p0.rect.width:.1f} x {p0.rect.height:.1f} pt "
      f"({p0.rect.width/72*25.4:.0f} x {p0.rect.height/72*25.4:.0f} mm)")
print(f"metadata   : title={doc.metadata.get('title')!r} "
      f"producer={doc.metadata.get('producer')!r}")

# ---------------------------------------------------------------- text scan
full = "\n".join(pg.get_text() for pg in doc)
print("\n" + "=" * 78)
print("CONTENT SANITY")
print("=" * 78)
for label, pat in (("unresolved refs '??'", r"\?\?"),
                   ("unresolved cites '[?]'", r"\[\?\]"),
                   ("leftover TODO", r"TODO"),
                   ("red placeholder text", r"author list|affiliation"),
                   ("'FILL IN' marker", r"FILL IN")):
    hits = re.findall(pat, full)
    print(f"  {label:26s} {len(hits)}")

print("\n--- first 30 non-empty lines of page 1 ---")
for ln in [l for l in doc[0].get_text().split("\n") if l.strip()][:30]:
    print("   ", ln[:100])

# ------------------------------------------------------------- column layout
print("\n" + "=" * 78)
print("PAGE LAYOUT")
print("=" * 78)
for pi, pg in enumerate(doc):
    blocks = [b for b in pg.get_text("blocks") if b[4].strip()]
    if not blocks:
        print(f"  p{pi+1}: (no text)")
        continue
    xs0 = min(b[0] for b in blocks)
    xs1 = max(b[2] for b in blocks)
    ys0 = min(b[1] for b in blocks)
    ys1 = max(b[3] for b in blocks)
    # detect a two-column layout by looking for a vertical gap in the middle
    mid = pg.rect.width / 2
    left = [b for b in blocks if b[2] <= mid + 6]
    right = [b for b in blocks if b[0] >= mid - 6]
    full_w = [b for b in blocks if b[0] < mid - 6 and b[2] > mid + 6]
    kind = "2-col" if (right and len(right) > 2) else "1-col"
    print(f"  p{pi+1}: {kind:5s} blocks={len(blocks):3d} "
          f"x[{xs0:6.1f},{xs1:6.1f}] y[{ys0:6.1f},{ys1:6.1f}] "
          f"left={len(left):3d} right={len(right):3d} spanning={len(full_w)}")

# ------------------------------------------------------- overflow detection
print("\n" + "=" * 78)
print("BLOCKS OUTSIDE THE PRINTABLE AREA  (Overfull hbox / vbox proxy)")
print("=" * 78)
# revtex PRE single column text width ~ 246 pt (3.4 in); margins ~ 54 pt
TOL = 2.0
worst = []
for pi, pg in enumerate(doc):
    W, H = pg.rect.width, pg.rect.height
    for b in pg.get_text("blocks"):
        if not b[4].strip():
            continue
        x0, y0, x1, y1, txt = b[0], b[1], b[2], b[3], b[4]
        over_r = x1 - (W - 40)
        over_l = 40 - x0
        over_b = y1 - (H - 30)
        if over_r > TOL or over_l > TOL or over_b > TOL:
            worst.append((max(over_r, over_l, over_b), pi + 1, x0, x1, y1,
                          txt.strip().replace("\n", " ")[:70]))
worst.sort(reverse=True)
if not worst:
    print("  none")
else:
    for w, pi, x0, x1, y1, txt in worst[:20]:
        print(f"  p{pi:3d} overflow={w:6.1f}pt  x=[{x0:6.1f},{x1:6.1f}] y1={y1:6.1f}"
              f"  {txt}")

# ------------------------------------------------------------ image / table
print("\n" + "=" * 78)
print("IMAGES PER PAGE")
print("=" * 78)
for pi, pg in enumerate(doc):
    imgs = pg.get_images(full=True)
    if imgs:
        rects = []
        for im in imgs:
            try:
                rects.extend(pg.get_image_rects(im[0]))
            except Exception:
                pass
        info = ", ".join(f"{r.width:.0f}x{r.height:.0f}pt@({r.x0:.0f},{r.y0:.0f})"
                         for r in rects)
        print(f"  p{pi+1}: {len(imgs)} image(s)  {info}")

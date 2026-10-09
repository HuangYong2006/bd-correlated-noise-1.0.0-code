r"""
Span-level overfull detection.

An Overfull \hbox shows up in the PDF as a text span whose bounding box extends
past its column edge.  RevTeX PRE in two-column mode has:
    text area   x in [54, 562] pt   (page width 612)
    left column  x in [54, 300]
    right column x in [316, 562]
so we flag:
    * right-column text     with x1 > 562 + tol
    * left-column text      with x1 > 300 + tol  (runs into the gutter)
    * anything              with x0 < 54 - tol
    * anything              with x1 > 562 + tol
and separately list the genuinely spanning material (figures/tables), which is
allowed and must be distinguished from overflow.

Usage: python code/overfull.py manuscript/overleaf_bd_paper.pdf
"""

from __future__ import annotations

import sys

import fitz

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

path = sys.argv[1] if len(sys.argv) > 1 else "manuscript/overleaf_bd_paper.pdf"
doc = fitz.open(path)

LEFT, MID, RIGHT = 54.0, 308.0, 562.0
TOL = 2.0

print(f"file: {path}   pages: {len(doc)}")
print(f"assumed text area x=[{LEFT:.0f},{RIGHT:.0f}], column split at {MID:.0f}")
print("=" * 100)

n_over = 0
per_page = {}
for pi, pg in enumerate(doc):
    d = pg.get_text("dict")
    for blk in d.get("blocks", []):
        if blk.get("type") != 0:
            continue
        for line in blk.get("lines", []):
            spans = line.get("spans", [])
            if not spans:
                continue
            x0 = min(s["bbox"][0] for s in spans)
            x1 = max(s["bbox"][2] for s in spans)
            txt = "".join(s["text"] for s in spans).strip()
            if not txt:
                continue
            # a line legitimately spanning both columns (title, abstract header,
            # figure*/table* caption) starts left of MID and ends right of it
            spans_both = x0 < MID - 6 and x1 > MID + 6
            over_right = x1 - RIGHT
            over_left = LEFT - x0
            into_gutter = (not spans_both) and x0 < MID - 6 and x1 > MID + 2
            if over_right > TOL or over_left > TOL or into_gutter:
                n_over += 1
                why = []
                if over_right > TOL:
                    why.append(f"past right margin by {over_right:.1f}pt")
                if over_left > TOL:
                    why.append(f"past left margin by {over_left:.1f}pt")
                if into_gutter:
                    why.append(f"runs into gutter by {x1-MID:.1f}pt")
                per_page.setdefault(pi + 1, []).append(
                    (max(over_right, over_left), x0, x1, txt[:88], "; ".join(why)))

if not per_page:
    print("No overfull lines detected at span level.")
for pi in sorted(per_page):
    print(f"\n--- page {pi} ---")
    for w, x0, x1, txt, why in sorted(per_page[pi], reverse=True):
        print(f"  x=[{x0:6.1f},{x1:6.1f}] {why}")
        print(f"      {txt}")

print("\n" + "=" * 100)
print(f"TOTAL overfull lines: {n_over}")
print("=" * 100)

# ---- what legitimately spans both columns --------------------------------
print("\nGenuinely spanning lines (figures / table* / title), for reference:")
for pi, pg in enumerate(doc):
    d = pg.get_text("dict")
    for blk in d.get("blocks", []):
        if blk.get("type") != 0:
            continue
        for line in blk.get("lines", []):
            spans = line.get("spans", [])
            if not spans:
                continue
            x0 = min(s["bbox"][0] for s in spans)
            x1 = max(s["bbox"][2] for s in spans)
            txt = "".join(s["text"] for s in spans).strip()
            if txt and x0 < MID - 6 and x1 > MID + 6:
                print(f"  p{pi+1:3d} x=[{x0:6.1f},{x1:6.1f}]  {txt[:80]}")

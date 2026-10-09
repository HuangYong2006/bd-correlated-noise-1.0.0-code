r"""
Shorten the header of the transposed alpha table.

The column headers were "FFGN (published)", "power-law filter" and
"exact fGn filter".  Each is 16 characters and at \footnotesize that is about
47 pt, so the three headers alone need 141 pt of a 246 pt column and the table
estimated at 91%.  Short labels bring it to about 45%, and the full names are
already given in the surrounding text.

Usage: python code/shorten_header.py
"""

from __future__ import annotations

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
p = os.path.join(ROOT, "manuscript", "main.tex")
src = open(p, encoding="utf-8").read()
n0 = len(src)

old = r"""$\theta$ & FFGN (published) & power-law filter & exact fGn filter\\"""
new = r"""$\theta$ & FFGN & power-law & exact fGn\\"""

c = src.count(old)
if c != 1:
    raise SystemExit(f"anchor count {c} (need 1)")
src = src.replace(old, new)

# make sure the surrounding text spells out the three generators
old2 = r"""realisations,"""
new2 = r"""realisations; the three columns are the FFGN generator, the power-law
spectral filter and the exact-fGn spectral filter,"""
c2 = src.count(old2)
if c2 == 1:
    src = src.replace(old2, new2)
    print("  lead-in sentence extended")
else:
    print(f"  (lead-in anchor appears {c2} times, left untouched)")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"header shortened; main.tex {n0} -> {len(src)} chars")

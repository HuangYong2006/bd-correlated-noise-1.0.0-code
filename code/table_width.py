r"""
Estimate the rendered width of every tabular in the manuscript, so that a table
that will overflow its column or the page can be caught WITHOUT compiling.

Method
  * parse each \begin{tabular}{<spec>} ... \end{tabular}
  * for l/c/r columns, take the maximum "visual length" of the cell text in that
    column over all rows (LaTeX commands stripped, braces removed, a few
    commands mapped to their rendered form, e.g. \alpha -> a, \theta -> q)
  * for p{...} / m{...} / b{...} columns, use the declared width and resolve
    \columnwidth / \textwidth numerically
  * width = sum(column widths) + 2*tabcolsep*ncols   (tabcolsep = 6 pt, revtex)
  * compare against the container: 246 pt for a single column, 508 pt for a
    full-width table* in revtex4-2 two-column mode on US Letter

The character-width constant is calibrated against a known case: the original
Appendix C table (plain tabular{ll}) measured 569.8 pt in the real compiled PDF.

Usage: python code/table_width.py manuscript/main.tex
"""

from __future__ import annotations

import os
import re
import sys

PT_PER_CHAR = 4.6          # calibrated (see module docstring)
TABCOLSEP = 6.0            # revtex default
# Target journal is Physical Review E, prepared with REVTeX 4.2 in the
# two-column layout: a single column is 246 pt and the full text width is 508 pt.
# The small inline reference tables sit inside ONE column; only table* floats get
# the full width.
COL_W = 246.0              # single column, REVTeX two-column
FULL_W = 508.0             # textwidth of the two-column layout

GLYPH = {
    r"\alpha": "a", r"\beta": "b", r"\gamma": "g", r"\theta": "q",
    r"\tau": "t", r"\eta": "h", r"\kappa": "k", r"\chi": "c", r"\nu": "v",
    r"\ell": "l", r"\Delta": "D", r"\Gamma": "G", r"\Omega": "W",
    r"\pm": "+-", r"\dots": "...", r"\ldots": "...", r"\times": "x",
    r"\sim": "~", r"\propto": "~", r"\to": "->", r"\ge": ">=",
    r"\leq": "<=", r"\pm0": "+-0", r"\alpha_s": "a_s",
}


def visual_len(s: str) -> int:
    """Rough rendered length of a LaTeX cell."""
    s = s.strip()
    for k in (r"\multicolumn", r"\cline", r"\hline", r"\centering"):
        s = s.replace(k, " ")
    # inline math: keep the content, drop the delimiters
    s = re.sub(r"\$([^$]*)\$", r"\1", s)
    for k, v in sorted(GLYPH.items(), key=lambda kv: -len(kv[0])):
        s = s.replace(k, v)
    s = re.sub(r"\\[A-Za-z]+", "x", s)      # any remaining command -> 1 glyph
    s = re.sub(r"[{}]", "", s)
    s = re.sub(r"\s+", " ", s)
    return len(s)


def parse_colspec(spec: str, ncols_hint: int):
    """Return a list of ('auto', None) or ('fixed', points) per column."""
    widths, i = [], 0
    while i < len(spec):
        ch = spec[i]
        if ch in "lcrX":
            widths.append(("auto", None)); i += 1
        elif ch in "pmb":
            m = re.match(r"[pmb]\{([^}]*)\}", spec[i:])
            if m:
                arg = m.group(1)
                if r"\columnwidth" in arg:
                    val = float(arg.replace(r"\columnwidth", "").strip("* ") or 1) * COL_W
                elif r"\textwidth" in arg:
                    val = float(arg.replace(r"\textwidth", "").strip("* ") or 1) * FULL_W
                else:
                    try:
                        val = float(arg)
                    except ValueError:
                        val = COL_W
                widths.append(("fixed", val))
                i += m.end()
            else:
                i += 1
        elif ch in "@!|><":
            m = re.match(r"[@!><]\{[^}]*\}", spec[i:])
            i += m.end() if m else 1
        elif ch == "*":
            i += 1
        else:
            i += 1
    if not widths:
        widths = [("auto", None)] * max(ncols_hint, 1)
    return widths


def rows_of(body: str):
    body = re.sub(r"\\hline|\\cline\{[^}]*\}", "", body)
    rows = []
    for raw in re.split(r"\\\\(?:\s*\[[^\]]*\])?", body):
        if not raw.strip():
            continue
        rows.append([c for c in raw.split("&")])
    return rows


def estimate(body: str, spec: str):
    cols = parse_colspec(spec, 1)
    n = len(cols)
    rows = rows_of(body)
    # multicolumn spans: widen the first covered column
    widths = [0.0] * n
    for r in rows:
        ci = 0
        for c in r:
            m = re.match(r"\s*\\multicolumn\{(\d+)\}", c)
            span = int(m.group(1)) if m else 1
            if ci >= n:
                break
            need = visual_len(c) * PT_PER_CHAR
            if span == 1:
                widths[ci] = max(widths[ci], need)
            else:
                span_need = need / span + (span - 1) * 2 * TABCOLSEP / span
                for k in range(span):
                    if ci + k < n:
                        widths[ci + k] = max(widths[ci + k], span_need)
            ci += span
    total = 0.0
    for k, (kind, val) in enumerate(cols):
        total += val if kind == "fixed" else widths[k]
    total += 2 * TABCOLSEP * n
    return total, n, [f"{w:.0f}" for w in widths]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "manuscript/main.tex"
    src = open(path, encoding="utf-8").read()

    print(f"estimator calibration constant: {PT_PER_CHAR} pt/char, "
          f"tabcolsep {TABCOLSEP} pt")
    print(f"column width {COL_W} pt, full width {FULL_W} pt\n")

    # ---- self test against the known 569.8 pt case ------------------------
    selftest_body = r"""
FFGN track (published generator) & $L=64\dots512$, $T=15L^{3/2}$, up to $256$ runs\\
spectral track (corrected) & $L=32\dots256$ ($512$ for selected $\theta$), $T=15L^{3/2}$, $192$ runs ($96$ at $L=512$)\\
Isolation Forest & $200$--$400$ trees, FPR $\alpha=0.02$ calibrated on held-out $\theta=0$ data, $3$--$16$ train/calibration resplits\\
"""
    est, n, _ = estimate(selftest_body, "ll")
    print(f"SELF TEST: the original Appendix C table (tabular{{ll}})")
    print(f"  estimated {est:6.1f} pt, measured in the real PDF 569.8 pt, "
          f"error {100*(est-569.8)/569.8:+.0f}%\n")

    # ---- every tabular in the file ---------------------------------------
    print("=" * 96)
    print(f"{'line':>6} {'container':>10} {'cols':>5} {'font':>9} {'sep':>5} "
          f"{'est width':>10} {'container':>10} {'use%':>6} {'fits?':>7}  label")
    print("-" * 96)
    bad, tight = 0, 0
    for m in re.finditer(r"\\begin\{tabular\}(?:\{([^}]*)\})?", src):
        spec = m.group(1) or ""
        start = m.end()
        endm = src.find(r"\end{tabular}", start)
        if endm < 0:
            continue
        body = src[start:endm]
        line = src[:m.start()].count("\n") + 1

        # font / tabcolsep overrides just before the tabular
        pre = src[max(0, m.start() - 600):m.start()]
        scale, fname = 1.0, "normal"
        if r"\compacttab" in pre.split("}")[-1] or r"\compacttab" in pre[-200:]:
            scale, fname = 0.80, "footnote"
        elif r"\footnotesize" in pre[-200:]:
            scale, fname = 0.80, "footnote"
        elif r"\small" in pre[-200:]:
            scale, fname = 0.90, "small"
        sep = TABCOLSEP
        sm = re.search(r"\\setlength\{\\tabcolsep\}\{([\d.]+)pt\}", pre[-200:])
        if sm:
            sep = float(sm.group(1))
        elif fname == "footnote" and r"\compacttab" in pre[-200:]:
            sep = 4.0

        before = src[max(0, m.start() - 4000):m.start()]
        if r"\begin{table*}" in before.split(r"\end{table*}")[-1]:
            container, limit = "table*", FULL_W
        elif r"\begin{table}" in before.split(r"\end{table}")[-1]:
            container, limit = "table", COL_W
        elif r"\begin{center}" in before.split(r"\end{center}")[-1]:
            container, limit = "centered", COL_W
        else:
            container, limit = "?", COL_W

        est_raw, n, _ = estimate(body, spec)
        # rescale: replace the auto-column part by scale, keep fixed columns
        cols = parse_colspec(spec, n)
        fixed = sum(v for k, v in cols if k == "fixed")
        auto = est_raw - fixed - 2 * TABCOLSEP * n
        est = fixed + auto * scale + 2 * sep * n
        pct = 100 * est / limit
        fits = est <= limit
        if not fits:
            bad += 1
        elif pct > 85:
            tight += 1
        rows = rows_of(body)
        label = rows[0][0][:26] if rows else ""
        print(f"{line:6d} {container:>10} {n:5d} {fname:>9} {sep:5.1f} "
              f"{est:10.1f} {limit:10.0f} {pct:5.0f}% "
              f"{'yes' if fits else 'NO':>7}  {label}")

    print("-" * 96)
    print(f"tables that do not fit : {bad}")
    print(f"tables above 85% of the container (risk margin): {tight}")
    print("\nEstimates are +/-20%.  Anything above 85% should be shortened.")


if __name__ == "__main__":
    main()

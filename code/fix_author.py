r"""
Final source cleanup, then rebuild of the Overleaf package.

main.tex still had the author/affiliation block wrapped in \TODO, which renders
red.  replace with plain placeholders.  The \TODO macro is kept defined (it is
useful for the author) and the header comment is updated to say that no \TODO
remains in the body.
"""

from __future__ import annotations

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
p = os.path.join(ROOT, "manuscript", "main.tex")
src = open(p, encoding="utf-8").read()
n0 = len(src)


def rep(old, new, label):
    global src
    c = src.count(old)
    if c != 1:
        raise SystemExit(f"anchor count {c} (need 1) for {label}")
    src = src.replace(old, new)
    print(f"  {label}")


rep(r"""%  NOTE ON STATUS: text marked with \TODO{...} is not yet final and is
%  rendered in red.  Every number that appears OUTSIDE a \TODO is backed by
%  a completed computation in the workspace (see README.md).""",
    r"""%  NOTE ON STATUS: no \TODO marker remains in the body.  Every number quoted
%  in the Results, Discussion and Conclusion is backed by a completed
%  computation in the workspace (see README.md).  \TODO is still defined and
%  renders in red should the author wish to mark further passages.""",
    "status comment updated")

rep(r"""\newcommand{\TODO}[1]{\textcolor{red}{[TODO: #1]}}""",
    r"""% red marker, available for the author's own use; currently unused
\newcommand{\TODO}[1]{\textcolor{red}{[TODO: #1]}}""",
    "TODO macro annotated")

rep(r"""\author{\TODO{author list}}
\affiliation{\TODO{affiliation}}""",
    r"""% ============ FILL IN ============
\author{First A. Author}
\affiliation{Department, Institution, City, Country}
\author{Second B. Author}
\affiliation{Department, Institution, City, Country}
% ==================================""",
    "author/affiliation placeholders")

open(p, "w", encoding="utf-8", newline="\n").write(src)
print(f"\n{n0} -> {len(src)} chars")
left = [f"L{i}" for i, ln in enumerate(src.split("\n"), 1)
        if "\\TODO{" in ln and not ln.strip().startswith("%")]
print(f"\\TODO{{...}} in the body: {left if left else 'none'}")

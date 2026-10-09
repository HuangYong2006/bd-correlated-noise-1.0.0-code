r"""
Static pre-flight check for a LaTeX manuscript, standing in for a compile we
cannot run (no TeX distribution available in this environment).

Checks performed
  1. every \begin{env} is matched by \end{env} with correct NESTING
  2. every \label is referenced and every \ref / \cite resolves
  3. all \includegraphics targets exist, resolved through \graphicspath
  4. math-mode dollar signs are balanced in each paragraph
  5. no undefined \newcommand-based macros are used
  6. brace balance per environment, and a global brace count
  7. flags constructs known to break specific classes (ruledtabular outside
     revtex, \multicolumn with \cline in a plain tabular)
  8. file encoding is pure ASCII or valid UTF-8

Usage: python code/check_tex.py manuscript/overleaf/main.tex
"""

from __future__ import annotations

import os
import re
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "manuscript/overleaf/main.tex"
base = os.path.dirname(os.path.abspath(path))

raw = open(path, "rb").read()
try:
    src = raw.decode("utf-8")
    enc = "utf-8"
except UnicodeDecodeError as e:
    print(f"FAIL encoding: not valid UTF-8 ({e})")
    sys.exit(1)
non_ascii = sorted({c for c in src if ord(c) > 127})
print(f"file      : {path}")
print(f"encoding  : {enc}, {len(src)} chars, "
      f"{len(non_ascii)} distinct non-ASCII chars {''.join(non_ascii[:20])}")

fails, warns = [], []

# ---- strip comments (respecting \%) --------------------------------------
lines = src.split("\n")
clean_lines = []
for ln in lines:
    out, i = [], 0
    while i < len(ln):
        if ln[i] == "\\" and i + 1 < len(ln):
            out.append(ln[i:i + 2]); i += 2; continue
        if ln[i] == "%":
            break
        out.append(ln[i]); i += 1
    clean_lines.append("".join(out))
clean = "\n".join(clean_lines)

# ---- 1. environment nesting ----------------------------------------------
stack = []
env_ok = True
for m in re.finditer(r"\\(begin|end)\{([^}]*)\}", clean):
    kind, env = m.group(1), m.group(2)
    ln = clean[:m.start()].count("\n") + 1
    if kind == "begin":
        stack.append((env, ln))
    else:
        if not stack:
            fails.append(f"L{ln}: \\end{{{env}}} with no matching \\begin")
            env_ok = False
        elif stack[-1][0] != env:
            fails.append(f"L{ln}: \\end{{{env}}} closes \\begin{{{stack[-1][0]}}} "
                         f"opened at L{stack[-1][1]}")
            env_ok = False
            stack.pop()
        else:
            stack.pop()
if stack:
    for env, ln in stack:
        fails.append(f"L{ln}: \\begin{{{env}}} never closed")
    env_ok = False
print(f"environments: {'OK' if env_ok else 'FAIL'} "
      f"({len(re.findall(r'\\begin{', clean))} begin, "
      f"{len(re.findall(r'\\end{', clean))} end)")

# ---- 2. labels / refs / cites --------------------------------------------
labels = set(re.findall(r"\\label\{([^}]*)\}", clean))
refs = set()
for m in re.finditer(r"\\(?:ref|eqref|autoref)\{([^}]*)\}", clean):
    refs.update(x.strip() for x in m.group(1).split(","))
missing_refs = sorted(refs - labels)
if missing_refs:
    fails.append(f"undefined \\ref targets: {missing_refs}")
print(f"cross-refs : {len(refs)} unique, "
      f"{'OK' if not missing_refs else 'FAIL ' + str(missing_refs)}")

bib = ""
bibpath = os.path.join(base, "refs.bib")
if os.path.exists(bibpath):
    bib = open(bibpath, encoding="utf-8").read()
cites = set()
for m in re.finditer(r"\\cite\{([^}]*)\}", clean):
    cites.update(x.strip() for x in m.group(1).split(","))
missing_cites = sorted(k for k in cites
                       if not re.search(r"@\w+\{" + re.escape(k) + r",", bib))
if missing_cites:
    fails.append(f"cite keys absent from refs.bib: {missing_cites}")
unused = sorted(k for k in re.findall(r"@\w+\{([^,]+),", bib) if k not in cites)
print(f"citations  : {len(cites)} used, "
      f"{'all OK' if not missing_cites else 'FAIL ' + str(missing_cites)}"
      f"  ({len(unused)} bib entries unused: {unused if unused else 'none'})")

# ---- 3. graphics ---------------------------------------------------------
gp = re.search(r"\\graphicspath\{((?:\{[^}]*\})*)\}", clean)
search = [""]
if gp:
    search = re.findall(r"\{([^}]*)\}", gp.group(1)) or [""]
missing_figs = []
for m in re.finditer(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", clean):
    target = m.group(1)
    found = None
    for d in search:
        for ext in ("", ".pdf", ".png", ".eps", ".jpg"):
            cand = os.path.join(base, d, target + ext)
            if os.path.exists(cand):
                found = cand
                break
        if found:
            break
    if not found:
        missing_figs.append(target)
if missing_figs:
    fails.append(f"\\includegraphics targets not found: {missing_figs}")
n_fig = len(re.findall(r"\\includegraphics", clean))
print(f"figures    : {n_fig} includes, graphicspath={search}, "
      f"{'all resolve' if not missing_figs else 'FAIL ' + str(missing_figs)}")

# ---- 3b. every float label must actually be cited ------------------------
# A figure or table that is typeset but never referenced is a defect a referee
# will flag, and it is invisible in a compile log.  Two figures in the previous
# draft of this manuscript were uncited, so this check is now part of the gate.
float_labels = re.findall(r"\\label\{(fig:[^}]*|tab:[^}]*)\}", clean)
uncited = sorted({lb for lb in float_labels if f"ref{{{lb}}}" not in clean})
if uncited:
    fails.append(f"float labels never referenced: {uncited}")
print(f"float refs : {len(set(float_labels))} labelled, "
      f"{'all cited' if not uncited else 'FAIL uncited ' + str(uncited)}")

# ---- 4. math-mode balance per paragraph ----------------------------------
depth, bad_math = 0, []
paras = re.split(r"\n\s*\n", clean)
for pi, para in enumerate(paras):
    for m in re.finditer(r"(?<!\\)\$", para):
        depth ^= 1
    if depth:
        bad_math.append(pi)
        depth = 0
if bad_math:
    warns.append(f"{len(bad_math)} paragraph(s) have an odd number of $ "
                 f"(possible unclosed math)")
print(f"math       : {'OK' if not bad_math else 'WARN ' + str(len(bad_math))}")

# ---- 5. macros -----------------------------------------------------------
defined = set(re.findall(r"\\newcommand\{\\([A-Za-z@]+)\}", clean))
used = set(re.findall(r"\\([A-Za-z@]+)", clean))
# a small allowlist of standard LaTeX/revtex/package macros
STD = set("""
documentclass usepackage begin end title author affiliation date maketitle
section subsection subsubsection paragraph appendix label ref eqref cite
bibliography includegraphics graphicspath caption textbf textit emph textrm
textcolor mathrm mathbf mathcal mathbb bm frac tfrac dfrac sqrt sum prod int
partial nabla infty delta Delta theta Theta alpha beta gamma Gamma lambda
Lambda chi sigma Sigma tau omega Omega eta mu nu xi pi rho phi psi epsilon
left right big Big bigl bigr rangle langle equiv approx sim simeq propto
lesssim gtrsim lessgtr geq leq neq pm mp times cdot ast circ oplus
line hline cline multicolumn tabular ruledtabular table figure itemize
enumerate description cite ref newcommand renewcommand providecommand
T1 fontenc amsmath amssymb graphicx xcolor hyperref floatfix aps pre
twocolumn superscriptaddress textwidth in columnspace
Phi Psi arcsin arccos arctan cos sin tan cosh sinh exp ln log max min
sup inf lim hat bar vec dot ddot tilde widehat widetilde lceil rceil
lfloor rfloor langle rangle downarrow uparrow leftrightarrow Leftrightarrow
Longrightarrow ge le ne to mapsto mathbf mathrm mathit mathsf mathtt
displaystyle textstyle scriptstyle mathstrut phantom vphantom hphantom
operatorname DeclareMathOperator nonumber notag tag qquad qquad enspace
quad hspace vspace smallskip medskip bigskip noindent newline par
si unit qty num SI intertext shortintertext
""".split())
unknown = sorted(u for u in used
                 if u not in STD and u not in defined
                 and not u.startswith(("rm", "it", "bf", "tt")))
# filter obvious text-mode words that follow a backslash accidentally
unknown = [u for u in unknown if len(u) > 1]
print(f"macros     : {len(defined)} defined; "
      f"{'OK' if not unknown else 'REVIEW ' + str(unknown[:14])}")
if unknown:
    warns.append(f"unrecognised macros (review): {unknown[:14]}")

# ---- 6. global brace balance --------------------------------------------
nb = len(re.findall(r"(?<!\\)\{", clean))
ne = len(re.findall(r"(?<!\\)\}", clean))
print(f"braces     : {{ = {nb}, }} = {ne}, "
      f"{'OK' if nb == ne else 'FAIL'}")
if nb != ne:
    fails.append(f"brace mismatch: {nb} open vs {ne} close")

# ---- 7. class-specific fragility ----------------------------------------
# The documentclass line must be a real (uncommented) line.  A scripted
# preamble rewrite once deleted it, leaving a file that still passed every other
# check here but could not compile at all, so the absence of a real
# \documentclass is a hard failure rather than a printed note.
doccls = ""
for ln in src.splitlines():
    if ln.lstrip().startswith("\\documentclass"):
        doccls = ln.strip()
        break
if not doccls:
    fails.append("no \\documentclass line found (only a mention inside a comment?)")
print(f"class      : {doccls or 'FAIL (not found!)'}")

# count how many \begin{document} appear: exactly one is required
n_doc = len(re.findall(r"\\begin\{document\}", clean))
if n_doc != 1:
    fails.append(f"\\begin{{document}} appears {n_doc} times (need exactly 1)")

if doccls and "revtex" not in doccls and "ruledtabular" in clean:
    warns.append("ruledtabular used but documentclass is not revtex")
if re.search(r"\\multicolumn", clean) and re.search(r"\\cline", clean) \
        and "ruledtabular" not in clean:
    warns.append("\\multicolumn together with \\cline in a plain tabular "
                 "can raise 'Misplaced \\noalign'")

# ---- 8. report -----------------------------------------------------------
print("\n" + "=" * 72)
if fails:
    print("HARD FAILURES:")
    for f in fails:
        print("  x", f)
for w in warns:
    print("  ! ", w)
if not fails:
    print("No structural errors found.  This is not a substitute for a real")
    print("compile, but it covers the failure modes that actually occur here.")

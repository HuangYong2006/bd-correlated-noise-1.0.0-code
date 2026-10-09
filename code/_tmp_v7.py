"""Cross-check the v7 PDF and the source for the remaining checklist items."""
import re
import collections
import fitz

PDF = r"C:\Users\MSI\Downloads\overleaf_bd_paper__7_.pdf"
doc = fitz.open(PDF)
txt = "\n".join(p.get_text() for p in doc)

print("=== P1-5: the glue strings the checklist lists ===")
for s in ["at\u03b8", "withL", "the\u03b8", "notwhite", "correlationamplitude",
          "exactpmf", "long-rangephase", "fixed-\u03b8scale", "crossoverlength"]:
    print(f"   {s!r:24} {'PRESENT' if s in txt else 'absent'}")
print("   spaced forms actually present in the PDF:")
for a, b in [("at", "\u03b8"), ("with", "L"), ("the", "\u03b8"),
             ("exact", "pmf"), ("crossover", "length"),
             ("correlation", "amplitude"), ("not", "white")]:
    n = len(re.findall(rf"{a}\s+{re.escape(b)}", txt))
    print(f"      '{a} {b}': {n}")

print("\n=== P1-8: reference list as rendered (entries 14-22) ===")
i = txt.find("REFERENCES")
if i < 0:
    i = txt.find("References")
tail = txt[i:]
parts = re.split(r"\n\s*\[(\d+)\]", "\n" + tail)
for k in range(1, len(parts), 2):
    n = int(parts[k])
    if 14 <= n <= 22:
        print(f"   [{n}] {' '.join(parts[k+1].split())[:230]}")

print("\n=== P1-6: source notation consistency ===")
src = open("manuscript/main.tex", encoding="utf-8").read()
body = "\n".join(l for l in src.splitlines() if not l.strip().startswith("%"))
for pat, label in [(r"\\mathrm\{nom\}", "\\mathrm{nom}"),
                   (r"\\rm nom", "\\rm nom"),
                   (r"\\mathrm\{th\}", "\\mathrm{th}"),
                   (r"\\rm th\b", "\\rm th"),
                   (r"\\mathrm\{eff\}", "\\mathrm{eff}"),
                   (r"\\rm eff", "\\rm eff"),
                   (r"\\mathrm\{flat\}", "\\mathrm{flat}"),
                   (r"\\rm flat", "\\rm flat"),
                   (r"\\mathrm\{sat\}", "\\mathrm{sat}"),
                   (r"\\rm sat", "\\rm sat")]:
    print(f"   {label:18} {len(re.findall(pat, body))}")

print("\n=== P1-6: L-range notation forms ===")
forms = collections.Counter(re.findall(r"L\s*=\s*32[^$]{0,14}", body))
for f, n in forms.most_common():
    print(f"   {f!r:26} x{n}")
print("   \\dots / -- usage around the range:")
print("     'L=32\\dots512' :", body.count(r"L=32\dots512"))
print("     '$L=32$--$512$' :", body.count(r"$L=32$--$512$"))
print("     '$L = 32$--$512$':", body.count(r"$L = 32$--$512$"))

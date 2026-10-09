r"""
Inspect a local Zotero library directly, without any plugin.

Zotero keeps everything in an ordinary SQLite database plus a storage tree:

    <data dir>/zotero.sqlite
    <data dir>/storage/<attachment item key>/<filename.pdf>

so anything with filesystem access can read the whole library -- metadata and
PDF full text -- using only the Python standard library plus a PDF extractor.

IMPORTANT implementation note: use a SEPARATE cursor for the per-item lookups.
Reusing the cursor that is iterating the outer query silently truncates the
outer result set to a single row; that bug produced the first, incomplete
listing of this library.

Usage
  python code/zotero_probe.py --read-pdfs 3
  python code/zotero_probe.py --list-collection 凝聚态 --limit 60
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

ap = argparse.ArgumentParser()
ap.add_argument("--datadir", default=os.path.join(os.path.expanduser("~"), "Zotero"))
ap.add_argument("--read-pdfs", type=int, default=5)
ap.add_argument("--list-collection", default=None)
ap.add_argument("--limit", type=int, default=40)
args = ap.parse_args()

DD = args.datadir
DB = os.path.join(DD, "zotero.sqlite")
STORAGE = os.path.join(DD, "storage")

print(f"data directory : {DD}")
if not os.path.exists(DB):
    raise SystemExit(f"database MISSING at {DB}")
print(f"database       : {DB}  ({os.path.getsize(DB)/1e6:.1f} MB)")

con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
cur = con.cursor()
sub = con.cursor()          # dedicated cursor for per-item lookups

for key in ("lastCompatibleVersion", "lastVersion"):
    try:
        v = sub.execute("SELECT value FROM settings WHERE key=?", (key,)).fetchone()
        if v:
            print(f"{key:15s}: {v[0]}")
    except sqlite3.Error:
        pass


def field(item_id, fname):
    r = sub.execute(
        "SELECT idv.value FROM itemData id "
        "JOIN fields f ON id.fieldID = f.fieldID "
        "JOIN itemDataValues idv ON id.valueID = idv.valueID "
        "WHERE id.itemID=? AND f.fieldName=?", (item_id, fname)).fetchone()
    return r[0] if r else ""


def first_author(item_id):
    r = sub.execute(
        "SELECT c.lastName FROM itemCreators ic JOIN creators c "
        "ON ic.creatorID = c.creatorID WHERE ic.itemID=? "
        "ORDER BY ic.orderIndex LIMIT 1", (item_id,)).fetchone()
    return r[0] if r else ""


print("\n" + "=" * 78)
print("LIBRARY CONTENTS")
print("=" * 78)
for name, n in cur.execute(
        "SELECT it.typeName, COUNT(*) FROM items i "
        "JOIN itemTypes it ON i.itemTypeID = it.itemTypeID "
        "GROUP BY it.typeName ORDER BY COUNT(*) DESC"):
    print(f"  {name:22s} {n:6d}")

REGULAR = ("journalArticle", "conferencePaper", "book", "bookSection",
           "thesis", "report", "preprint", "manuscript", "webpage")
ph = ",".join("?" * len(REGULAR))
n_regular = sub.execute(
    f"SELECT COUNT(*) FROM items i JOIN itemTypes t ON i.itemTypeID=t.itemTypeID "
    f"WHERE t.typeName IN ({ph})", REGULAR).fetchone()[0]
print(f"\n  -> bibliographic items: {n_regular}")

print("\n" + "=" * 78)
print("COLLECTIONS")
print("=" * 78)
for name, n in cur.execute(
        "SELECT c.collectionName, COUNT(ci.itemID) FROM collections c "
        "LEFT JOIN collectionItems ci ON c.collectionID = ci.collectionID "
        "GROUP BY c.collectionID ORDER BY COUNT(ci.itemID) DESC"):
    print(f"  {n:6d}  {name}")

if args.list_collection:
    print("\n" + "=" * 78)
    print(f"ITEMS IN COLLECTION matching '{args.list_collection}'")
    print("=" * 78)
    rows = cur.execute(
        f"SELECT i.itemID FROM items i "
        f"JOIN collectionItems ci ON i.itemID = ci.itemID "
        f"JOIN collections c ON ci.collectionID = c.collectionID "
        f"JOIN itemTypes t ON i.itemTypeID = t.itemTypeID "
        f"WHERE c.collectionName LIKE ? AND t.typeName IN ({ph}) "
        f"ORDER BY i.dateAdded DESC LIMIT ?",
        (f"%{args.list_collection}%", *REGULAR, args.limit)).fetchall()
else:
    print("\n" + "=" * 78)
    print(f"MOST RECENT {args.limit} ITEMS")
    print("=" * 78)
    rows = cur.execute(
        f"SELECT i.itemID FROM items i "
        f"JOIN itemTypes t ON i.itemTypeID = t.itemTypeID "
        f"WHERE t.typeName IN ({ph}) ORDER BY i.dateAdded DESC LIMIT ?",
        (*REGULAR, args.limit)).fetchall()

for (item_id,) in rows:
    t = field(item_id, "title") or "(no title)"
    y = (field(item_id, "date") or field(item_id, "year"))[:4]
    doi = field(item_id, "DOI")
    print(f"  {first_author(item_id):16s} {y:4s}  {t[:70]}")
    if doi:
        print(f"                           doi:{doi}")

print("\n" + "=" * 78)
print("ATTACHED PDFs")
print("=" * 78)
pdfs, missing = [], 0
for att_id, att_key, path, parent in cur.execute(
        "SELECT a.itemID, i.key, a.path, a.parentItemID FROM itemAttachments a "
        "JOIN items i ON a.itemID = i.itemID "
        "WHERE a.contentType = 'application/pdf'"):
    if not path:
        continue
    fn = path.split(":", 1)[1] if path.startswith("storage:") else path
    full = os.path.join(STORAGE, att_key, os.path.basename(fn))
    if os.path.exists(full):
        pdfs.append((att_id, parent, full))
    else:
        missing += 1
print(f"  linked PDF attachments : {len(pdfs) + missing}")
print(f"  present on disk        : {len(pdfs)}")
print(f"  missing on disk        : {missing}")
if pdfs:
    print(f"  total size             : "
          f"{sum(os.path.getsize(p) for _, _, p in pdfs)/1e6:.0f} MB")

if args.read_pdfs and pdfs:
    print("\n" + "=" * 78)
    print(f"TEXT EXTRACTION TEST on {min(args.read_pdfs, len(pdfs))} PDFs")
    print("=" * 78)
    try:
        import fitz
    except ImportError:
        print("  PyMuPDF unavailable: pip install pymupdf")
        raise SystemExit(0)

    ok = 0
    for att_id, parent, path in pdfs[:args.read_pdfs]:
        title = field(parent, "title") if parent else ""
        try:
            doc = fitz.open(path)
            n_chars = sum(len(doc[i].get_text()) for i in range(min(3, len(doc))))
            head = " ".join((doc[0].get_text() if len(doc) else "").split())[:92]
            good = n_chars > 400
            ok += good
            print(f"\n  {os.path.basename(path)[:56]}")
            print(f"    pages={len(doc):4d}  chars(first3)={n_chars:7d}  "
                  f"{'TEXT OK' if good else 'LITTLE TEXT (scan?)'}")
            if title:
                print(f"    title : {title[:66]}")
            print(f"    p1    : {head}")
            doc.close()
        except Exception as e:
            print(f"\n  {os.path.basename(path)}\n    FAILED {type(e).__name__}: {e}")
    print(f"\n  -> {ok}/{min(args.read_pdfs, len(pdfs))} readable as text")

con.close()
print("\nDONE")

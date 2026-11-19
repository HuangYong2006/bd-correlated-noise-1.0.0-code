"""Render complete PDFs and record mechanical publication-artifact checks."""
import hashlib
import json
from pathlib import Path
import re

import fitz
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT/'docs/research_refocus_20261117'
BUILD = DOCS/'pdf_output'


def main():
    preview = BUILD/'preview_final'
    preview.mkdir(exist_ok=True)
    records = []
    for name, count in [('main', 5), ('supplemental', 4)]:
        path = BUILD/f'{name}.pdf'
        document = fitz.open(path)
        text = '\n'.join(p.get_text() for p in document)
        assert '??' not in text and '\ufffd' not in text
        log = (BUILD/f'{name}.log').read_text(errors='replace')
        assert not re.search(r'Overfull|undefined references|undefined citations|Missing character', log, re.I)
        prefix = '' if name == 'main' else 'S'
        for i in range(1, count+1):
            assert f'FIG. {prefix}{i}.' in text
        thumbs, fonts = [], {}
        for i, page in enumerate(document):
            for f in page.get_fonts(full=True):
                assert f[2] != 'Type3'
                if f[0]:
                    fonts[f[3]] = len(document.extract_font(f[0])[3]) > 0
            for word in page.get_text('words'):
                assert word[0] >= -1 and word[1] >= -1
                assert word[2] <= page.rect.width+1 and word[3] <= page.rect.height+1
            pix = page.get_pixmap(matrix=fitz.Matrix(1.7, 1.7), alpha=False)
            pix.save(preview/f'{name}_page{i+1}.png')
            picture = Image.frombytes('RGB', [pix.width, pix.height], pix.samples)
            picture.thumbnail((510, 660))
            tile = Image.new('RGB', (530, 690), '#eeeeee')
            tile.paste(picture, ((530-picture.width)//2, 20))
            ImageDraw.Draw(tile).text((12, 4), f'{name} page {i+1}', fill='black')
            thumbs.append(tile)
        assert all(fonts.values()), fonts
        contact = Image.new('RGB', (1590, 690*((len(thumbs)+2)//3)), 'white')
        for i, tile in enumerate(thumbs):
            contact.paste(tile, ((i%3)*530, (i//3)*690))
        contact.save(preview/f'{name}_contact.png')
        records.append(dict(file=path.relative_to(ROOT).as_posix(), pages=len(document),
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                            metadata=document.metadata, embedded_fonts=fonts,
                            figure_count=count, unresolved_references=False, overfull_boxes=False))
    (DOCS/'pdf_qa.json').write_text(json.dumps(records, indent=2)+'\n', encoding='utf-8')
    print('PDF checks passed; inspect rendered contact sheets and individual pages before delivery.')


if __name__ == '__main__':
    main()

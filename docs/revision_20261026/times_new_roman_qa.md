# Figure font revision

- Scope: the four figures included by the current PRE manuscript. Python plotting workflow retained. Existing quantitative panels, data, run selections, statistics, colours and output dimensions retained.
- All plot text uses Times New Roman, including regular, bold, italic and mathematical labels. Embedded PDF font names were checked for all four figures and again recursively inside the compiled manuscript. No Arial or DejaVu fonts were found in the compiled manuscript.
- Regenerated PDF, SVG and 600 dpi PNG files. The source manifest was identical before and after regeneration, including input SHA-256 values and representative-run indices.
- Figure 1 left margin increased from 0.075 to 0.08 to accommodate the font change and prevent its first vertical label from clipping. All panels were visually inspected.
- Exported-PDF text audits passed for all four figures. Minimum text sizes were 5.04, 5.25, 5.25 and 5.04 pt respectively.
- Source validator FONT-FAMILY finding is overridden by the user's explicit Times New Roman request. The validator's default expects sans-serif publication styling. Actual embedded fonts satisfy the requested family. No validator rules were changed.
- Source validator warnings were reviewed. Actual rendered glyphs pass the size floor. The width is 178 mm, not the validator's incorrectly inferred 4521.2 mm. SVG/PDF are vector deliverables and PNG is the preview, so no TIFF was added for this font-only revision.
- The existing manuscript compiled successfully with the editor compiler and its exported PDF was regenerated. The PDF has eight pages and no unresolved-reference markers.

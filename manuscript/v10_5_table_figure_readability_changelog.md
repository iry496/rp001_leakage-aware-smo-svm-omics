# v10.5 — Table & Figure Readability Cleanup — Changelog

**Source:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_4_final_copyedit_consistency.docx`
**Output:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_5_table_figure_readability.docx`
**Date:** 2026-06-24
**Scope:** Layout/formatting only. No analysis; no result values changed; no new results; no new references. A whole-document text comparison confirms the body and table text are byte-identical to v10.4 (only formatting properties changed).

## Decisions (confirmed before editing)
- **Tables:** compact only — no consolidation / no master matrix (lowest risk; preserves all granular tables and values).
- **Dashboard figure:** rotate to its own landscape page.

## Changes made

### 1. Table compaction (Tasks 1, 2, 3) — applied to all 13 tables
- **Root cause fixed:** the document default line spacing is double (line = 480), which was also being applied *inside table cells*, making every row unnecessarily tall and pushing tables across extra pages with large whitespace gaps.
- Set **single line spacing** inside all table cells (409 cell paragraphs), with 0 pt before / 2 pt after. This compacts every table, especially the priority tables (2, 5, 6, 7, 10, 11, 12).
- **Table 12 (Reproducible Omics Evidence Audit)** previously spanned three pages; it now fits in about two with no broken-looking rows.
- Header-row repeat (`tblHeader`) and "don't split row across pages" (`cantSplit`) were already present on the tables and are retained — so headers still repeat on continuation pages and no row breaks mid-cell.
- The 11-column Selector Top-K Sweep table (Table 6) already fit the page width; it remains intact and is now less cramped vertically.
- Overall document length: **49 → 47 pages.**
- **No consolidation performed** (per decision). Tables 5, 6, 10, 11, and 12 remain separate; no master performance matrix was created. No values moved or merged.

### 2. Figure rotation (Tasks 5, 6) — the 6-panel Evidence Audit dashboard
- **Numbering note:** the dense 6-panel dashboard is **Figure 8** in the manuscript (caption: "Figure 8. Reproducible Omics Evidence Audit dashboard…"). The task referred to it as "Figure 7"; Figure 7 is actually the clean three-group transportability bar chart, which is already readable and was left unchanged. The change below applies to Figure 8.
- Placed Figure 8 on **its own landscape page** via two section breaks (a portrait section ending just before the figure, a landscape section ending after the caption, then portrait resumes). Only that single page is landscape; all other pages remain portrait.
- **Enlarged Figure 8** from 6.50 in to **9.00 in** wide (height scaled proportionally to 4.95 in; aspect ratio preserved, source image is high-resolution 2942×1617). All six panels are noticeably larger and more legible.
- The figure caption travels with the figure on the landscape page.
- Page numbering was kept continuous across the new sections (removed the page-number restart so numbering still runs 1…47 without resetting after the landscape page). Footer/page number retained on the landscape page.

## Verification
- DOCX validates and renders (245 paragraphs, 13 tables; rendered to PDF without error).
- **All body and table text identical to v10.4** (formatting-only change) — confirmed by full-document text diff.
- **All numeric tokens identical to v10.4** — no result value changed.
- Figures **1–8** and Tables **1–12** remain sequential; Figure 8 keeps its number (not renumbered, since it was rotated rather than split).
- Page-orientation check: only the Figure 8 page is landscape; all surrounding pages are portrait.
- GSE25065 / GSE41998 framing, leakage-vs-transportability distinction, and all prior consistency properties are unaffected (no text changed).

# Data Availability Statement

All datasets analyzed in this study are publicly available from the NCBI Gene Expression Omnibus (GEO; https://www.ncbi.nlm.nih.gov/geo/). No new primary data were generated. Expression matrices were obtained as author-processed GEO Series Matrix files via the `GEOquery` package (`GSEMatrix = TRUE`, `getGPL = FALSE`); raw CEL files were not reprocessed.

## Datasets used

| Accession | Platform | Role in this study | Source publication |
| --- | --- | --- | --- |
| GSE25055 | Affymetrix HG-U133A (GPL96) | Discovery / internal nested cross-validation | Hatzis et al., 2011 (JAMA; GSE25066 SuperSeries) |
| GSE25065 | Affymetrix HG-U133A (GPL96) | Same-platform, same-study-family external validation | Hatzis et al., 2011 (JAMA; GSE25066 SuperSeries) |
| GSE41998 | Affymetrix HG-U133A 2.0 (GPL571) | Predeclared cross-platform transportability sensitivity | Horak et al., 2013 (Clin Cancer Res) — **citation pending author approval (see citation verification report)** |

GSE25055 and GSE25065 are the non-overlapping discovery and validation subseries of the GSE25066 SuperSeries.

## Datasets documented but not used

GSE20194 and GSE20271 were documented for audit transparency but were **not** included in the reported analysis, owing to MDACC-lineage patient-overlap risk that would require sample-level de-duplication before use.

## Code and reproducibility artifacts

The analysis code, configuration, random seeds, fold assignments, selected feature lists, committed performance tables, and external-validation outputs are available in the project repository. The repository is publicly available at https://github.com/iry496/leakage-aware-smo-svm-omics. A frozen DOI-citable archive will be deposited before final journal submission.

> No Zenodo DOI is asserted because none has been minted yet. Once the archive is deposited, insert the DOI here and in `CITATION.cff`, the manuscript Availability statement, and the release notes. **[AUTHOR/RELEASE ACTION: mint Zenodo DOI and backfill.]**

## Ethics

This study used only publicly available, de-identified, previously published human gene-expression data obtained from GEO. No new human-subjects data were collected and no identifiable patient information was accessed; institutional ethics approval was therefore not required for this secondary computational analysis. Ethical approvals for the original cohorts are described in their respective source publications.

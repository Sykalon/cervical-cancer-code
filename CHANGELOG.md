# Preparation changes (unreleased)

- Removed machine-specific `setwd()` paths; run each analysis from its own input directory.
- Renamed the immune infiltration script to `IOBR.R` for a portable repository filename.
- WGCNA: commented pasted console output so the script parses; replaced an unimported pipe with base function-call syntax; omitted the interactive graphics window.
- RF: corrected the malformed `ggsave()` call and explicitly saved the first survival plot; explicitly printed the confusion plot to its PDF device; corrected the top-10 comment to match the actual top-5 code.
- NMF: explicitly loaded survival and survminer; saved the rank survey as `rank_survey.pdf` to avoid overwriting it with the coefficient map.
- Added documentation, source checksums and ignore patterns for local data/results.
- Scientific parameters, rank selection, sample-count assumptions, method choices and time units were not silently changed. Unresolved items are documented in REVIEW_REQUIRED.md.

## Author update

- Hui Ye supplied revised scripts: WGCNA now uses `n_samples` for correlation p-values; the undefined `cnt` print was removed from RF.
- Author metadata recorded as Hui Ye (GitHub: Sykalon).

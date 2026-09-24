# cervical-cancer-code

Deconstructing Cervical Cancer Heterogeneity: Squamous Cell Carcinoma Subtypes and CCM2 Vulnerability in Adenocarcinoma

Source-code archive containing four R analysis scripts: WGCNA, immune infiltration estimation with IOBR, random forest classification, and NMF clustering.

**Status:** Not a validated reproduction release. See [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md) for unresolved settings and runtime issues. This package contains code only, not patient data or trained models.

## Files and required inputs

Use a separate working directory for each analysis: several scripts use the same input/output filenames for different data.

| Script | Input files | Expected structure from supplied code |
| --- | --- | --- |
| WGCNA.R | data.csv, trait.csv | Expression: features in rows, samples in columns, first column row identifiers. Traits: samples in rows, numeric traits in columns, first column sample identifiers. Sample alignment must be checked. |
| IOBR.R | ydata.csv | Features in rows, samples in columns, first column feature identifiers. Identifier type and preprocessing need author confirmation. |
| RF.R | scale.csv, data.csv, PFS.csv | Training: samples in rows, numeric predictors and Group label. External data: matching predictors and name identifier. PFS: name, Days, Status. Time units and event coding need confirmation. |
| NMF.R | sd.csv, 总生存率.csv | Nonnegative feature-by-sample matrix with first column feature identifiers. Survival: ID, days, Status. Time units and event coding need confirmation. |

## Dependencies and execution

Packages referenced by these scripts: WGCNA, stringr, IOBR, randomForest, caret, dplyr, survival, survminer, reshape2, ggprism, ggplot2, scales, NMF, doMPI. Installation requirements and study package versions have not yet been established. An environment lockfile cannot be reconstructed from scripts alone.

After resolving the review items and preparing inputs, run a script from its own working directory, for example:

```r
setwd("path/to/rf-working-directory")
source("path/to/repository/RF.R", encoding = "UTF-8")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

Each script reads relative input filenames and writes outputs in the working directory. Existing output files may be overwritten. Scripts clear the R workspace; use a fresh R session. There is no integrated pipeline or bundled dataset.

## Validation and provenance

Only syntax checks are performed during publication preparation. Analyses have not been rerun and scientific outputs have not been verified. SOURCE_MANIFEST.json records checksums of the supplied source files. Changes are listed in CHANGELOG.md. Original supplied files are retained locally outside this publication package.

## Citation and license

Author: Hui Ye. GitHub account: Sykalon. Repository: https://github.com/Sykalon/cervical-cancer-code. Licensed under the MIT License; copyright (c) 2026 Hui Ye. CITATION.cff records the software citation. Zenodo archival details will be added when available.


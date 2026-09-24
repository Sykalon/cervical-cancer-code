# Validation and known limitations

This release archives analysis scripts. Syntax parsing passed; the study analyses have not been rerun during release preparation. Input data and the original package environment are not included. Scientific settings were retained from the supplied scripts.

- WGCNA: verify trait/sample alignment and missing values. Module correlations default to Pearson. The secondary TOM uses power 4 whereas the main network uses an estimated power. Hub and network exports select different modules. The final export references annotation columns not created earlier and an undefined `design` object; that section requires adaptation before execution.
- Random forest: external-column intersection may reduce the five selected predictors. Preprocessing of scale.csv is not documented here. Cross-validation is commented out; the active code reports accuracy. Confusion plot axis labels are reversed relative to the predicted/reference table and percentages are normalized within predicted classes. Verify survival time units (Days versus month labels) and event coding. SHAP analysis is not included in these supplied scripts.
- NMF: review rank selection from the cophenetic differences. The script contains repeated fits, a method comparison and a final rank-4 fit that overwrites the model object. The main fit does not explicitly specify its algorithm. Verify survival time units and the doMPI environment. The TIFF device is closed without drawing a plot.
- IOBR: the script includes eight methods, including TIMER with a coad label. Confirm method/input compatibility, feature identifiers, normalization and the cancer-type setting before reuse.

These limitations are recorded to distinguish a source-code archive from a fully validated reproduction pipeline.

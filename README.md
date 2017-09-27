CamHack17
================

Single-cell RNA-Seq cell-type prediction
-----------------------------------------
R package to predict the cell types of new single-cell RNA-Seq samples based on past collection of annotated gold-standard single-cell RNA-Seq experiments:
- Supply a collection of annotated samples. 33 public annotated gold-standard scRNA-Seq datasets are already collected by Lin et al. [Using Neural Networks To Improve Single-Cell RNA-Seq Data Analysis](https://www.biorxiv.org/content/early/2017/04/23/129759) -- this collection can be used as a starting point.
- Predict new samples using a basic PCA-based model

Future direction:
- Incorporate more advanced methods such as deep learning
- Show prediction accuracy using e.g. MAP scores
- Cell-type ontology-aware predictions

Automatic P-value thresholding
-----------------------------------------
In Bioinformatics, we do multiple hypothesis testing a lot, e.g. in GWAS and differential expression analysis. This introduces many false positives and currently, we perform a multiple-testing correction and choose an arbitrary threshold to limit them. But, this manual thresholding is sub-optimal and your significant results may unnecessarily include too many false positives or exclude too many true positives. An R package will be implemented for an automatic and data-driven alternative to finding an optimal P-value threshold.

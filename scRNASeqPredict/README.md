scRNASeqPredict Guide
================
Kevin Dialdestoro
2017-09-27

Use several annotated gold standard scRNA-Seq datasets for training and testing the prediction model, as already collected by Lin et al. [Using Neural Networks To Improve Single-Cell RNA-Seq Data Analysis](https://www.biorxiv.org/content/early/2017/04/23/129759). Datasets are collected from [here](http://sb.cs.cmu.edu/scnn/).

Read training datasets and their metadata
-----------------------------------------

``` r
base.path = '/Users/kevindialdestoro/scRNASeqPredict'
# datasets: GEO accessions for spleen (x2), ESC, thymus (x2), kidney, lung
datasets.accessions = c('GSE75109', 'GSE75110',
                        'E-MTAB-2805',
                        'GSE60297', 'GSE74596',
                        'GSE66202',
                        'GSE66578')
datasets.paths = paste(base.path, 'data', paste0(datasets.accessions, '.txt'), sep='/')
metadata.path = paste(base.path, 'metadata.tsv', sep='/')

exprs = ReadDatasets(datasets.paths)
metadata = ReadMetadata(metadata.path)
table(metadata[colnames(exprs),'celltype']) # summary of cell types across the datasets
```

    ## 
    ##                     ESC_G1                    ESC_G2M 
    ##                         96                         96 
    ##                      ESC_S kidney_nephron_progenitors 
    ##                         96                         91 
    ##          lung_CD4+_T_cells         spleen_and_LN_Th17 
    ##                          6                        269 
    ##                thymus_mTEC                thymus_NKT0 
    ##                        174                         45 
    ##                thymus_NKT1                thymus_NKT2 
    ##                         90                         68

Build PCA model
---------------

``` r
pca.model = BuildPCAModel(exprs)
```

Test prediction on a dataset with spleen cells
----------------------------------------------

``` r
# test 1 using dataset GSE75111 of 151 spleen cells
query.dataset.accession = 'GSE75111'
query.dataset.path = paste(base.path, 'data', paste0(query.dataset.accession, '.txt'), sep='/')
query.expr = read.table(query.dataset.path, sep="\t", header=TRUE, check.names=FALSE, row.names=1)
predictions = PredictCellTypes(query.expr, pca.model, metadata)
table(predictions) # 100% accuracy
```

    ## predictions
    ## spleen_and_LN_Th17 
    ##                151

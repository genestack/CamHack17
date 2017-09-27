library(devtools)
library(flashpcaR)

ReadDatasets = function(...) {
    datasets.paths = c(...)
    expr.list = lapply(datasets.paths, read.table, sep="\t", header=TRUE, check.names=FALSE, row.names = NULL, stringsAsFactors = FALSE)

    # merge tables into one by common genes
    # transpose s.t. the samples correspond to rows
    expr = Reduce(function(...) merge(..., by='gene'), expr.list)
    rownames(expr) = expr[,'gene']
    expr = expr[,-1]
    expr = data.matrix(expr)
    return(expr)
}

ReadMetadata = function(metadata.path) {
    metadata = read.table(metadata.path, sep="\t", header=TRUE, check.names=FALSE, row.names = 1, stringsAsFactors = FALSE)
    return(metadata)
}

BuildPCAModel = function(expr) {
    # remove genes with all zeroes
    expr = expr[rowSums(expr) != 0, ]

    number.components = min(100, ncol(expr))
    pca.model = flashpca(expr, ndim=number.components, stand="sd")

    # PCA weights
    pcaWeights = pca.model$vectors
    rownames(pcaWeights) = rownames(expr)
    colnames(pcaWeights) = paste0("PC",1:number.components)

    # PCA values
    pcaValues = t(expr) %*% pca.model$vectors
    rownames(pcaValues) = colnames(expr)
    colnames(pcaValues) = paste0("PC",1:number.components)
    rownames(pcaValues) = paste0(rownames(pcaValues))

    return(list(weights = pcaWeights, values = pcaValues))
}

PredictCellTypes = function(query.expr, pca.model, metadata) {
    common.genes = row.names(pca.model$weights)
    query.expr = query.expr[common.genes, ]

    predictions = sapply(1:ncol(query.expr), function(i) {
        query.expr.single = query.expr[, i]
        testCellPCA = t(query.expr.single) %*% pca.model$weights
        distances = apply(pca.model$values, 1, function(x){dist(rbind(x,testCellPCA))})
        distances = sort(distances)
        closestCells = names(head(distances,10))
        closestCells = metadata[closestCells, ]
        return(names(sort(table(closestCells),decreasing=TRUE)[1]))
    })
    names(predictions) = colnames(query.expr)
    return(predictions)
}

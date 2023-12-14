#!/usr/bin/env R

# Author: Sean Maden
#
# Test music2Param using lute example data.
#
#
#

library(lute)

data <- lute:::getDeconvolutionExampleData_music2()

singleCellExpressionSet <- data$singleCellExpressionSet
singleCellExpressionSet[["experimentCondition"]] <- 
  c(rep("control", 100), rep("disease", 200))

singleCellExperiment <- SummarizedExperiment(scEset) |> se_to_sce()
names(assays(singleCellExperiment)) <- "counts"
assays(singleCellExperiment)[["counts"]] <- 
  exprs(assays(singleCellExperiment)[["counts"]])
colData(singleCellExperiment) <- pData(singleCellExpressionSet) |> DataFrame()
referenceExpression <- 
  referenceFromSingleCellExperiment(singleCellExperiment, "counts", "cellType")

bulkExpressionSet <- data$bulkExpressionSet
bulkExpressionSet$experimentCondition <- 
  c(rep("control", 6), rep("disease", 6))
bulkExpression <- exprsMatrix <- exprs(bulkExpressionSet)
scFilterSamples <- 
  !colnames(bulkExpression) %in% singleCellExpressionSet$SubjectName
bulkExpressionIndependent <- bulkExpression[,scFilterSamples]
cellScaleFactors <- rep(1, length(unique(singleCellExpressionSet$cellType)))
sceNew <- SummarizedExperiment(assays = exprsMatrix)

newParam <- music2Param(bulkExpression = bulkExpression, 
                        bulkExpressionIndependent = bulkExpressionIndependent, 
                        referenceExpression = matrix(), 
                        cellScaleFactors = cellScaleFactors, 
                        singleCellExpressionSet = singleCellExpressionSet,
                        singleCellExperiment = singleCellExperiment)

names(attributes(newParam))
object <- newParam
deconvolution(newParam)






data <- lute:::.get_decon_example_data_bisque()
z <- lute:::.get_z_from_sce(sce, "counts", "cellType")
y.eset <- data$y.eset
exprs.matrix <- exprs(y.eset)
y <- exprs.matrix
filter.sc.samples <- !colnames(y) %in% sc.eset$SubjectName
yi <- y[,filter.sc.samples]
s <- rep(1, length(unique(sc.eset$cellType)))

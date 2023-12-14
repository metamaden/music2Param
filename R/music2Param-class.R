#!/usr/bin/env R

# Author: Sean Maden

#' music2Param-class
#' 
#' Main constructor for class to manage mappings to the `MuSiC2` deconvolution 
#' algorithm implementations \code{MuSiC::music2_prop() and 
#' \code{MuSiC2::music2_prop()}}. For the latter you will need to have 
#' installed an older version of MuSiC (< v1.0.0) as of MuSiC2 v0.1.0.
#' 
#' @details Main constructor for class \linkS4class{music2Param}.
#' @rdname music2Param-class
#' @seealso 
#' \linkS4class{referencebasedParam-class},
#' \linkS4class{independentbulkParam}
#' 
#' @aliases 
#' music2Param-class, MuSiC2Param-class, Music2Param-class
#' 
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- music2Param(exampleDataList[["bulkExpression"]],
#' exampleDataList[["referenceExpression"]],
#' exampleDataList[["cellScaleFactors"]])
#' deconvolution(newParam)
#' 
#' @references 
#' 
#' Fan, Jiaxin. MuSiC2: cell type deconvolution for multi-condition bulk 
#' RNA-seq data. (2023) GitHub, R package version 0.1.0. URL: 
#' https://github.com/Jiaxin-Fan/MuSiC2
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#' 
#' Jiaxin Fan, Yafei Lyu, Qihuang Zhang, Xuran Wang, Mingyao Li, Rui Xiao, 
#' MuSiC2: cell-type deconvolution for multi-condition bulk RNA-seq data, 
#' Briefings in Bioinformatics, Volume 23, Issue 6, November 2022, bbac430, 
#' https://doi-org.proxy1.library.jhu.edu/10.1093/bib/bbac430
#'
setClass("music2Param", 
         contains="independentbulkParam", 
         slots=c(
           singleCellExpressionSet = "ExpressionSet", 
           singleCellExperiment = "SingleCellExperiment",
           assayName = "character", 
           batchVariable = "character",
           cellTypeVariable = "character", 
           cellSizeScaleFactors = "numeric", 
           cellTypeSubset = "character",
           conditionVariable = "character",
           controlLabel = "character",
           caseLabel = "character",
           methodType = "character")
         )

#' Make an object of class music2Param
#'
#' Main function to make a new object of class \linkS4class{music2Param}, with
#' defaults for required arguments.
#'
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param bulkExpressionIndependent  Bulk mixed signals matrix of independent 
#' samples, which should not overlap samples in y.
#' @param referenceExpression Signature matrix of cell type-specific signals. 
#' If not provided, can be computed from a provided ExpressionSet containing 
#' single-cell data.
#' @param cellScaleFactors  Cell size factor transformations of length equal to 
#' the K cell types to deconvolve.
#' @param SingleCellExpressionSet ExpressionSet of single-cell expression.
#' @param SingleCellExperiment SingleCellExperiment of single-cell expression.
#' @param assayName Expression data type (e.g. counts, logcounts, tpm, etc.).
#' @param batchVariable Name of variable identifying the batches in 
#' SingleCellExpressionSet pData/coldata.
#' @param cellTypeVariable Name of cell type labels variable in 
#' SingleCellExpressionSet pData/coldata.
#' @param conditionVariable Name of variable in y.eset and 
#' SingleCellExpressionSet containing condition labels.
#' @param controlLabel Label of control condition samples in condition variable.
#' @param caseLabel Label of case condition samples in condition variable.
#' @param methodType Name of method source library to call for music2_prop; 
#' either "MuSiC" or "MuSiC2".
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#'
#' @seealso \linkS4class{musicParam}
#'
#' @returns Object of class \linkS4class{music2Param}.
#' 
music2Param <- function(bulkExpression = NULL, bulkExpressionIndependent = NULL, 
                        referenceExpression = NULL, cellScaleFactors = NULL, 
                        singleCellExpressionSet = NULL, 
                        singleCellExperiment = NULL,
                        assayName = "counts", batchVariable = "SubjectName", 
                        cellTypeVariable = "cellType", 
                        conditionVariable = "experimentCondition", 
                        controlLabel = "control", caseLabel = "case",
                        methodType = "MuSiC2", returnInfo = FALSE) {
  new("music2Param", 
      bulkExpression = bulkExpression, 
      bulkExpressionIndependent = bulkExpressionIndependent, 
      referenceExpression = referenceExpression, 
      cellScaleFactors = cellScaleFactors, 
      singleCellExpressionSet = singleCellExpressionSet, 
      singleCellExperiment = singleCellExperiment,
      assayName = assayName, batchVariable = batchVariable, 
      cellTypeVariable = cellTypeVariable, 
      conditionVariable = conditionVariable, 
      controlLabel = controlLabel, caseLabel = caseLabel, 
      methodType = methodType, 
      returnInfo = returnInfo)
}

#' Deconvolution method for \linkS4class{music2Param}
#'
#' Main method to access the MuSiC2 deconvolution algorithm.
#'
#' @param object An object of class \linkS4class{music2Param}.
#'
#' @details Takes an object of class \linkS4class{music2Param} as input, 
#' returning a list or vector of predicted cell type proportions.
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#' 
#' @references 
#' 
#' Fan, Jiaxin. MuSiC2: MuSiC2: cell type deconvolution for multi-condition bulk 
#' RNA-seq data. (2023) GitHub, R package version 0.1.0. URL: 
#' https://github.com/Jiaxin-Fan/MuSiC2
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#' 
#' Jiaxin Fan, Yafei Lyu, Qihuang Zhang, Xuran Wang, Mingyao Li, Rui Xiao, 
#' MuSiC2: cell-type deconvolution for multi-condition bulk RNA-seq data, 
#' Briefings in Bioinformatics, Volume 23, Issue 6, November 2022, bbac430, 
#' https://doi-org.proxy1.library.jhu.edu/10.1093/bib/bbac430
#'
setMethod("deconvolution", signature(object = "music2Param"), function(object){
  require(Biobase)
  # load data
  lparam <- callNextMethod()
  object <- lparam[["object"]]
  
  # instantiate function objects
  bulkExpressionIndependent <- object[["bulkExpressionIndependent"]]
  bulkExpression <- object[["bulkExpression"]]
  singleCellExpressionSet <- object[["singleCellExpressionSet"]]
  singleCellExperiment <- object[["singleCellExperiment"]]
  cellTypeSubset <- object[["cellTypeSubset"]]
  batchVariable <- object[["batchVariable"]]
  conditionVariable <- object[["conditionVariable"]]
  
  
  cellTypeVariable <- object[["cellTypeVariable"]]
  cellScaleFactors <- data.frame(
    cellTypes = singleCellExperiment[[cellTypeVariable]] |> unique(),
    cellSizeScaleFactors = object[["cellScaleFactors"]]
  )
  uniqueTypes <- cellScaleFactors[,1]
  
  caseLabel <- object[["caseLabel"]]
  controlLabel <- object[["controlLabel"]]
  # get result according to method type
  sourceLibrary <- object[["methodType"]]
  if(sourceLibrary %in% c("music", "MuSiC", "Music", "MUSIC")){
  	message("Using the MuSiC implementation of music2_prop()...")
    require(MuSiC)
  	result <- MuSiC::music2_prop(bulk.control.mtx = bulkExpression, 
  	                             bulk.case.mtx = bulkExpressionIndependent, 
  	                             sc.sce = singleCellExperiment, 
  	                             clusters = cellTypeVariable, 
  	                             samples = batchVariable, 
  	                             cell_size = cellScaleFactors, 
  	                             select.ct = NULL)
  } else{
  	message("Using the MuSiC2 implementation of music2_prop()...")
    require(MuSiC2)
	  result <- MuSiC2::music2_prop(
	    bulk.eset = bulkExpressionSet, 
	    sc.eset = singleCellExpressionSet, 
      condition = conditionVariable, 
      control = controlLabel, 
      case = caseLabel, 
      clusters = cellTypeVariable,
      samples = batchVariable, 
      cell_size = cellScaleFactors, 
      select.ct = uniqueTypes)
  }
  # return results
  predictions <- matrix(
    result$bulk.props, ncol = ncol(z))
  predictions <- apply(predictions, 1, function(ri){ri/sum(ri)})
  colnames(predictions) <- colnames(z)
  rownames(predictions) <- colnames(y)
  returnList <- t(predictions)
  if(object[["returnInfo"]]){
    returnList <- list(
      predictions = predictions, 
      result.info = result, 
      metadata = list(lmd = lparam[["metadata"]])
    )
  }
  return(returnList)
})

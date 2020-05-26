

#' Function to summarise the results of multiple phyllosphere experiments
#' @param cdsList A list object containing compatible data sets created by "makeCDS". The names of each item in the list will be used to label the summary tables. The results in the data sets must contain the same rows in the same order, i.e.: this function is designed to summarise peturbations to a core experiment.
#' @details
#' None
#' @keywords None
#' @return A list containing two data frames: fcMatrix is the summary of log2FoldChange results, pvMatrix is the summary of adjusted p-value results (padj).
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

summariseResults <- function(cdsList){
    resultsList <- lapply(cdsList,function(x) x$results)

    fcMatrix <- lapply(resultsList,function(x) x[,"log2FoldChange",F])
    fcMatrix <- do.call(cbind,fcMatrix)
    colnames(fcMatrix) <- names(resultsList)
    fcMatrix[is.na(fcMatrix)] <- 0

    pvMatrix <- lapply(resultsList,function(x) x[,"padj",F])
    pvMatrix <- do.call(cbind,pvMatrix)
    colnames(pvMatrix) <- names(resultsList)
    pvMatrix[is.na(pvMatrix)] <- 1

    return(list(fcMatrix=fcMatrix,pvMatrix=pvMatrix))
}

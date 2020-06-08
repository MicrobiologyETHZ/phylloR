

#' Function to produce a pie chart of an individual community.
#'
#' @param counts A count table in which the rows are strains or OTUs and the columns are samples. The table should include only a single sample group.
#' @param strainTaxa A vector of the taxonomic labels of the strains, in the order they appear in the rows of the count table.
#' @param cols A vector of colours where the names() of the vector correspond to the taxonomic labels.
#' @param taxLabels An optional vector defining the taxonomic labels to plot, with the default being the unique values of strainTaxa.
#' @param sort.tax If true, the taxonomic labels will be sorted from smallest to largest abundance.
#' @param drop A logical indicating whether unseen taxa should be dropped from the plot, defaults to TRUE.
#' @param ... Further arguments to be passed to "pie".
#' @details
#' For each taxonomic label, the sum of all counts for the strains with that label across all samples is divided by the total sum of the count table.
#' @keywords phylloR
#' @return None.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

# Pie chart of a community at a provided taxonomic level, using totals across multiple samples
plotCommunityPie <- function(counts,strainTaxa,cols=NULL,taxLabels=NULL,sort.tax=FALSE,drop=TRUE,...){
    if(is.null(cols)){
        cols <- rainbow(32)
    }
    if(is.null(taxLabels)){
        taxLabels <- unique(strainTaxa)
    }

    strainTotals <- apply(counts,1,sum)
    taxTotals <- sapply(taxLabels,function(x) sum(strainTotals[strainTaxa==x]))

    if(sort.tax){
        taxTotals <- sort(taxTotals)
        taxLabels <- names(taxTotals)
    }
    if(drop){
        taxLabels <- taxLabels[taxTotals>0]
        taxTotals <- taxTotals[taxTotals>0]
    }

    pie(taxTotals,radius=0.4,col=cols[taxLabels],clockwise=T,...)
    legend(0.3,0.4,xjust=0,yjust=0,taxLabels,fill=cols[taxLabels])
}

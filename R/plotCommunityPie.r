#' Function to produce a pie chart of an individual community.
#'
#' @param counts A count table in which the rows are strains or OTUs and the columns are samples. The table should include only a single sample group.
#' @param level  The taxonomic level at which the data will be summarised
#' @param cols   A vector of colours where the names() of the vector correspond to the taxonomic labels.
#' @details
#' For each taxon at the chosen level, the sum of all counts for the strains in that taxon across all samples is divided by the total sum of the count table.
#' @keywords phylloR
#' @return None.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

# Pie chart of a community at a provided taxonomic level, using totals across multiple samples
plotCommunityPie <- function(counts,level="Phlass",cols=phlassColors){
    # Fetch the relevant taxonomy table and add a line for unknowns
    validLeafs <- rownames(counts)[rownames(counts)%in%rownames(leafTaxonomy)]
    relTaxonomy <- leafTaxonomy[validLeafs,]

    # Normalise and summarise table
    normCounts <- apply(counts,2,function(x) x/sum(x))
    validCounts <- normCounts[validLeafs,,drop=F]

    taxLevels <- levels(as.factor(relTaxonomy[,level]))
    taxCounts <- sapply(taxLevels,function(x) apply(validCounts[relTaxonomy[,level]==x,],2,sum,na.rm=T))
    finalCounts <- t(cbind(taxCounts,Unclassified=1-apply(taxCounts,1,sum)))

    pie(finalCounts,radius=0.4,col=c(cols[taxLevels],"gray"),clockwise=T,...)
    legend(0.3,0.4,xjust=0,yjust=0,taxLevels,fill=c(cols[taxLevels],"gray"))
}

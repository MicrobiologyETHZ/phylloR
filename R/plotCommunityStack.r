#' Function to produce a stacked bar chart of a phyllosphere community
#'
#' @param counts    A count table in which the rows are strains or OTUs and the columns are samples.
#' @param level     Taxonomic level at which to summarise the community
#' @param taxLevels Optional vector of taxonomic levels, in case you only want to plot some or reorder them
#' @param meta      Optional metadata that should be used to group samples
#' @param cols      A named vector of colors for the different taxonomic levels
#' @param mode      Either 'relative' or 'absolute' for calculating abundance
#' @details
#'
#' @keywords phylloR
#' @return None.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

# Stacked bar chart of a phyllosphere community
plotCommunityStack <- function(counts,level="Phlass",taxLevels=NULL,meta=NULL,cols=phlassColors,mode="relative"){
    if(!mode%in%c("relative", "absolute")){
        cat("Mode not valid, choose from 'relative' or 'absolute'\n")
        return()
    }

    # Fetch the relevant taxonomy table and add a line for unknowns
    validLeafs <- rownames(counts)[rownames(counts)%in%rownames(leafTaxonomy)]
    relTaxonomy <- leafTaxonomy[validLeafs,,drop=F]

    # Normalise and summarise table
    normCounts <- apply(counts,2,function(x) x/sum(x))
    validCounts <- normCounts[validLeafs,,drop=F]

    if(is.null(taxLevels)){
        taxLevels <- levels(as.factor(relTaxonomy[,level]))
    }
    taxCounts <- sapply(taxLevels,function(x) apply(validCounts[relTaxonomy[,level]==x,,drop=F],2,sum,na.rm=T))
    finalCounts <- t(cbind(taxCounts,Unclassified=1-apply(taxCounts,1,sum)))

    if(!is.null(meta)){    
        finalCounts <- sapply(unique(meta),function(x) apply(finalCounts[,meta==x],1,sum))
        if(mode=="relative"){
            finalCounts <- apply(finalCounts,2,function(x) x/sum(x))
        }
    }

    par(mar=c(10,4,4,2)+0.1,xpd=T)
    barplot(finalCounts,col=c(cols[taxLevels],"gray"),names.arg=unique(meta))
    legend("bottom",inset=c(0,-0.2),horiz=T,legend=taxLevels,fill=c(cols[taxLevels],"gray"))
}

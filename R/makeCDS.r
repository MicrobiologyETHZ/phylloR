#' Function to make a compatible data set (CDS) and calculate DESeq2 results.
#'
#' @param counts A count table in which the rows are strains or OTUs and the columns are samples.
#' @param meta A metadata table in which the rows correspond to the columns of the count table.
#' @param include A list that indicates which samples should be included in the analysis, see details.
#' @param exclude A list that indicates which samples should be excluded from the analysis, see details.
#' @param foi The factor of interest that determines how samples will be grouped for comparison.
#' @param ftc The factor to consider separately from the feature of interest, for instance to account for a batch effect.
#' @param title A title for the data set that will carry through to plots.
#' @param legend A vector with the names of the sample groups.
#' @details
#' As an alternative to giving both the *counts* and *meta* arguments, you may give a single argument that is a list containing items labelled *counts* and *meta*.
#' The lists for the include and exclude arguments should consist of vectors named to match the columns of the metadata being considered and contain the valid (or invalid) items.
#' For instance: include=list(initial_condition=c("Control","Test"),time=c("t1","t2")).
#' Samples are filtered first by the include list, then the exclude list. You can of course slice the data set in advance according to more complex criteria.
#' @keywords phylloR
#' @return A list containing the original dataset, the DESeq2-analysed dataset and the accompanying results.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

makeCDS <- function(counts,meta,include=NULL,exclude=NULL,foi,ftc=NULL,title=NULL,legend=NULL){
    if(missing(meta)){
        ds = counts
        counts = ds$counts
        meta = ds$meta
    }
    # Slice data
    if(!is.null(include)){
        for(i in 1:length(include)){
            counts = counts[,meta[,names(include)[i]]%in%include[[i]]]
            meta = meta[meta[,names(include)[i]]%in%include[[i]],]
        }
    }
    if(!is.null(exclude)){
        mask <- rep(TRUE,ncol(counts))
        for(i in 1:length(exclude)){
            mask = mask & meta[,names(exclude)[i]]%in%exclude[[i]]
        }
        counts = counts[,!mask]
        meta = meta[!mask,]
    }
    
    # Remove unused factor levels and factorise key metadata
    for(i in 1:ncol(meta)){
        if(is.factor(meta[,i])){
            meta[,i] <- droplevels(meta[,i])
        }
    }
    if(!is.factor(meta[,foi])){
        meta[,foi] <- as.factor(meta[,foi])
    }
    if(!is.null(ftc)){
        if(!is.factor(meta[,ftc])){
            meta[,ftc] <- as.factor(meta[,ftc])
        }
    }

    # Run DESeq
    if(!is.null(ftc)){
        dds <- DESeqDataSetFromMatrix(counts,meta,as.formula(paste("~",ftc,"+",foi)))
    }else{
        dds <- DESeqDataSetFromMatrix(counts,meta,as.formula(paste("~",foi)))
    }
    dds <- DESeq(dds)
    results <- as.data.frame(results(dds))

    # Make compatible data set
    cds <- list(counts=counts,meta=meta,foi=foi,ftc=ftc,title=title,legend=legend,dds=dds,results=results)
    return(cds)
}

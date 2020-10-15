#' Function to produce a volcano plot of DESeq2 results
#'
#' @param cds       A compatible data set produced by the makeCDS() function.
#' @param cutoff    A p-value cutoff to determine which strains will be included individually in the plot.
#' @param pcols     An optional vector of point colors
#' @param pcex      An optional numeric vector indicating the size of the points
#' @param pborder   An optional vector of border colors for the points
#' @param hili      Logical; whether to highlight significant points or not
#' @param labels    Logical; whether to label significant points or not
#' @param ranks     A vector of ranks for the points in the plot that will be labelled accordingly
#' @details
#' None.
#' @keywords phylloR
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

plotVolcano <- function(cds,cutoff=0.05,pcols=NULL,pcex=NULL,pborder=NULL,hili=TRUE,labels=TRUE,ranks=NULL){
    if(is.null(pcols)){
        pcols = leafTaxonomy[rownames(cds$results),]$Color
        pcols[is.na(pcols)] <- "#BEBEBE"
    }
    if(is.null(pcex)){
        pcex = 2+4*cds$results$baseMean/max(cds$results$baseMean)
    }
    if(is.null(ranks)){
        ranks = order(cds$results$baseMean)
    }

    ylimit <- 10^floor(log10(min(cds$results$padj[!is.na(cds$results$padj)])))
    xlimit <- ceiling(max(abs(cds$results$log2FoldChange[!is.na(cds$results$log2FoldChange)])))
    psig <- which(cds$results$padj<=cutoff)

    if(hili){
        alphas <- c("77", "FF")[1+(cds$results$padj<=cutoff)]
        alphas[is.na(alphas)] <- "77"
        pcols <- paste(pcols, alphas, sep="")
    }

    if(is.null(pborder)){
        plot(cds$results$log2FoldChange,cds$results$padj,xlim=c(-xlimit,xlimit),ylim=c(1,ylimit),log="y",xlab="Log2 Fold Change",ylab="Adjusted P-Value",pch=20,col=pcols,cex=pcex,panel.first=grid())
    }else{
        plot(cds$results$log2FoldChange,cds$results$padj,xlim=c(-xlimit,xlimit),ylim=c(1,ylimit),log="y",xlab="Log2 Fold Change",ylab="Adjusted P-Value",pch=21,col=pborder,bg=pcols,cex=pcex,panel.first=grid())
    }
    abline(h=cutoff,lty=2)
    if(labels & (length(psig)>0)){
        text(cds$results$log2FoldChange[psig],cds$results$padj[psig],rownames(cds$results)[psig],pos=4,offset=0.2)
    }
    if(length(psig)>0){
        text(cds$results$log2FoldChange[psig],cds$results$padj[psig],ranks[psig],adj=0.5,cex=0.5)
    }
}

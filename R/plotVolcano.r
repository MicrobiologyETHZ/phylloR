#' Function to produce a volcano plot of DESeq2 results
#'
#' @param cds       A compatible data set produced by the makeCDS() function.
#' @param cutoff    A p-value cutoff to determine which strains will be included individually in the plot.
#' @param hili       Logical; whether to highlight significant points or not
#' @param label     Logical; whether to label significant points or not
#' @details
#' None.
#' @keywords phylloR
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

plotVolcano <- function(cds,cutoff=0.05,hili=TRUE,label=TRUE){
    ylimit <- 10^floor(log10(min(cds$results$padj[!is.na(cds$results$padj)])))
    xlimit <- ceiling(max(abs(cds$results$log2FoldChange[!is.na(cds$results$log2FoldChange)])))
    psig <- which(cds$results$padj<=cutoff)

    if(hili){
        pcols <- 1+(cds$results$padj<=cutoff)
    }else{
        pcols <- rep(1,nrow(cds$results$padj))
    }

    plot(cds$results$log2FoldChange,cds$results$padj,xlim=c(-xlimit,xlimit),ylim=c(1,ylimit),log="y",xlab="Log2 Fold Change",ylab="Adjusted P-Value",pch=20,col=pcols,panel.first=grid())
    abline(h=cutoff,lty=2)
    if(label){
        text(cds$results$log2FoldChange[psig],cds$results$padj[psig],rownames(cds$results)[psig],pos=4,offset=0.2)
    }
}

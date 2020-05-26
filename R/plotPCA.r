

#' Function to calculate and plot a PCA comparing sample groups at the whole community level.
#'
#' @param cds A compatible data set produced by the makeCDS() function.
#' @param soi A vector containing the names of the strains of interest, matching row names of the CDS count table.
#' @param perm The number of permutations for the PERMANOVA, passed to the adonis() function.
#' @param cutoff A p-value cutoff to determine which strains will be included individually in the plot.
#' @param rowLabs A vector of names for the strains of interest, for plotting.
#' @param subtitle A subtitle for the plot.
#' @param cols A vector of 2 colours for the points.
#' @param showLegend A logical indicating whether or not to show a legend.
#' @param showArrows A logical indicating whether or not to show the top 3 or fewer significant individual responses
#' @param showTitle A logical indicating whether or not to show a plot title.
#' @details
#' None.
#' @keywords phylloR
#' @return The results of the PCA and the adonis() function.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

plotPCA <- function(cds,soi=NULL,perm=100,cutoff=0.05,rowLabs=NULL,subtitle=NULL,cols=1:2,showLegend=TRUE,showArrows=TRUE,showTitle=TRUE){
    if(is.null(soi)){
        soi = rownames(cds$counts)
    }
    if(is.null(rowLabs)){
        rowLabs = soi
    }
    if(is.null(subtitle)){
        title = cds$title
    }else{
        title = paste(cds$title,subtitle,sep="\n")
    }

    if(showLegend & showTitle){
        par(mar=0.1+c(7,4,4,9),xpd=T)
    }else if(showLegend){
        par(mar=0.1+c(7,4,1,9),xpd=T)
    }else if(showTitle){
        par(mar=0.1+c(7,4,4,1))
    }else{
        par(mar=0.1+c(7,4,1,1))
    }

    nct <- assay(varianceStabilizingTransformation(cds$dds,blind=F))
    nct <- t(nct)[,soi]
    pca <- prcomp(nct)
    if(!is.null(cds$ftc)){
        adn <- adonis(as.formula(paste("nct~",cds$foi)),strata=cds$meta[,cds$ftc],perm=perm,cds$meta,method="euclidean")
    }else{
        adn <- adonis(as.formula(paste("nct~",cds$foi)),perm=perm,cds$meta,method="euclidean")
    }

    if(!is.null(cds$ftc)){
        plot(pca$x[,1:2],
             col=cols[cds$meta[,cds$foi]],
             pch=14+as.numeric(cds$meta[,cds$ftc]),
             main=title,
             cex=1.5
            )
        title(sub=paste("Effect Size: ",formatC(100*adn$aov.tab$R2[1],digits=3),"%; P-value: ",formatC(adn$aov.tab$Pr[1],3),sep=""),line=5)
    }else{
        plot(pca$x[,1:2],
             col=cols[cds$meta[,cds$foi]],
             pch=16,
             main=title,
             cex=1.5
            )
        title(sub=paste("Effect Size: ",formatC(100*adn$aov.tab$R2[1],digits=3),"%; P-value: ",formatC(adn$aov.tab$Pr[1],3),sep=""),line=5)
    }

    if(showArrows){
        res <- cds$results[soi,]
        res <- cbind(res,label=rowLabs)
        res <- res[order(res$padj),]
        res <- res[res$padj<cutoff & !is.na(res$padj),,drop=F]
        if(nrow(res)>0){
            ntop <- min(nrow(res),3)
            top <- res[1:ntop,]
            arrows(rep(0,ntop),rep(0,ntop),pca$rotation[rownames(top),1]*10,pca$rotation[rownames(top),2]*10,length=0.1,col=c("blue","white","red")[sign(top[,2])+2])
            text(pca$rotation[rownames(top),1]*15,pca$rotation[rownames(top),2]*15,top$label,col=c("blue","white","red")[sign(top[,2])+2])
        }
    }

    if(showLegend){
        if(!is.null(cds$ftc)){
            pairs = expand.grid(cds$legend,levels(cds$meta[,cds$ftc]))
            legend("topleft",inset=c(1.01,0),legend=apply(pairs,1,paste,collapse=" "),pch=14+as.numeric(pairs[,2]),col=rep(cols,length(levels(cds$meta[,cds$ftc]))))
        }else{
            legend("topleft",inset=c(1.01,0),legend=cds$legend,pch=16,col=cols)
        }
    }
    return(list(pca=pca,stats=adn))
}

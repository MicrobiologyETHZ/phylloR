

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
#' @param showHull A logical indicating whether or not to highlight the hull of each set of points
#' @param showSidebars A logical indicating whether or not to add box-and-whisker plots per group beside the axes
#' @param showTitle A logical indicating whether or not to show a plot title.
#' @param calcSpread    A logical indicating whether or not to calculate the spread of each group of points
#' @details
#' The spread of points is calculated as the average distance to the mean location of those points.
#' @keywords phylloR
#' @return The results of the PCA and the adonis() function, and optionally the spread of points in each group.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

plotPCA <- function(cds,soi=NULL,perm=100,cutoff=0.05,rowLabs=NULL,subtitle=NULL,cols=1:2,showLegend=TRUE,showArrows=TRUE,showHull=FALSE,showSidebars=FALSE,showTitle=TRUE,calcSpread=TRUE){
    if(is.null(soi)){
        soi = rownames(cds$counts)
    }
    if(is.null(rowLabs)){
        rowLabs = soi
    }
    if(is.null(subtitle)){
        plotTitle = cds$title
    }else{
        plotTitle = paste(cds$title,subtitle,sep="\n")
    }

    # Margin adjustment
    if(showSidebars){
        nSidebars = length(cols)
    }else{
        nSidebars = 0
    }
    par(mar=0.1+c(7,4,1+(3*showTitle)+(2*nSidebars),1+(8*showLegend)+(2*nSidebars)),xpd=(showLegend | showSidebars))

    nct <- assay(varianceStabilizingTransformation(cds$dds,blind=F))
    nct <- t(nct)[,soi]
    pca <- prcomp(nct)
    if(!is.null(cds$ftc)){
        adn <- adonis(as.formula(paste("nct~",cds$foi)),strata=cds$meta[,cds$ftc],perm=perm,cds$meta,method="euclidean")
    }else{
        adn <- adonis(as.formula(paste("nct~",cds$foi)),perm=perm,cds$meta,method="euclidean")
    }

    # Calculate variance explained
    varxp <- 100*pca$sdev^2/sum(pca$sdev^2)

    if(!is.null(cds$ftc)){
        plot(pca$x[,1:2],
             col=cols[cds$meta[,cds$foi]],
             pch=14+as.numeric(cds$meta[,cds$ftc]),
             cex=1.5,
             xlab = paste("PC1 (",round(varxp[1],3),")",sep=""),
             ylab = paste("PC2 (",round(varxp[2],3),")",sep="")
            )
        title(sub=paste("Effect Size: ",formatC(100*adn$aov.tab$R2[1],digits=3),"%; P-value: ",formatC(adn$aov.tab$Pr[1],3),sep=""),line=5)
    }else{
        plot(pca$x[,1:2],
             col=cols[cds$meta[,cds$foi]],
             pch=16,
             cex=1.5,
             xlab = paste("PC1 (",round(varxp[1],3),")",sep=""),
             ylab = paste("PC2 (",round(varxp[2],3),")",sep="")
            )
        title(sub=paste("Effect Size: ",formatC(100*adn$aov.tab$R2[1],digits=3),"%; P-value: ",formatC(adn$aov.tab$Pr[1],3),sep=""),line=5)
    }
    if(showTitle){
        title(main=plotTitle)
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

    if(showHull){
        sets = split(as.data.frame(pca$x[,1:2]),cds$meta[,cds$foi])
        for(i in 1:length(sets)){
            set = sets[[i]]
            hcoords = set[chull(set),]
            polygon(hcoords,col=adjustcolor(cols[i],alpha.f=0.2))
        }
    }

    if(showSidebars){
        # PC1
        sets = split(pca$x[,1],cds$meta[,cds$foi])
        scale = (max(pca$x[,2])-min(pca$x[,2]))/20
        for(i in 1:length(sets)){
            set = sets[[i]]
            box = boxplot(set,plot=F)
            lines(c(box$stats[1],box$stats[5]),rep(max(pca$x[,2])+((i+0.4)*scale),2),col=cols[i])
            rect(box$stats[2],max(pca$x[,2])+(i*scale),box$stats[4],max(pca$x[,2])+((i+0.8)*scale),border=cols[i],col="white")
            lines(rep(box$stats[3],2),c(max(pca$x[,2])+(i*scale),max(pca$x[,2])+((i+0.8)*scale)),col=cols[i])
            points(box$out,rep(max(pca$x[,2])+((i+0.4)*scale),length(box$out)),pch=19,cex=0.5,col=cols[i])
        }
        # PC2
        sets = split(pca$x[,2],cds$meta[,cds$foi])
        scale = (max(pca$x[,1])-min(pca$x[,1]))/20
        for(i in 1:length(sets)){
            set = sets[[i]]
            box = boxplot(set,plot=F)
            lines(rep(max(pca$x[,1])+((i+0.4)*scale),2),c(box$stats[1],box$stats[5]),col=cols[i])
            rect(max(pca$x[,1])+(i*scale),box$stats[4],max(pca$x[,1])+((i+0.8)*scale),box$stats[2],border=cols[i],col="white")
            lines(c(max(pca$x[,1])+(i*scale),max(pca$x[,1])+((i+0.8)*scale)),rep(box$stats[3],2),col=cols[i])
            points(rep(max(pca$x[,1])+((i+0.4)*scale),length(box$out)),box$out,pch=19,cex=0.5,col=cols[i])
        }

    }

    if(showLegend){
        scale = (max(pca$x[,1])-min(pca$x[,1]))/20
        if(!is.null(cds$ftc)){
            pairs = expand.grid(cds$legend,levels(cds$meta[,cds$ftc]))
            legend(max(pca$x[,1])+((1+nSidebars)*scale),par('usr')[4],legend=apply(pairs,1,paste,collapse=" "),pch=14+as.numeric(pairs[,2]),col=rep(cols,length(levels(cds$meta[,cds$ftc]))))
        }else{
            legend(max(pca$x[,1])+((1+nSidebars)*scale),par('usr')[4],legend=cds$legend,pch=16,col=cols)
        }
    }

    if(calcSpread){
        sets = split(as.data.frame(pca$x), cds$meta[, cds$foi])
        spreads <- c()
        for(i in 1:length(sets)){
            centroid = apply(sets[[i]], 2, mean)
            distances = apply(sets[[i]], 1, function(x) sqrt(sum((x-centroid)^2)))
            spreads[levels(cds$meta[, cds$foi])[i]] <- mean(distances)
        }
    }


    return(list(pca=pca, stats=adn, spreads=spreads))
}

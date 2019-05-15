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

#' Function to plot the changes in the community at a strain level.
#'
#' @param cds A compatible data set produced by the makeCDS() function.
#' @param soi A vector containing the names of the strains of interest, matching row names of the CDS count table.
#' @param cutoff A p-value cutoff to determine which strains will be included individually in the plot.
#' @param rowLabs A vector of names for the strains of interest, for plotting.
#' @param subtitle A subtitle for the plot.
#' @param cols A vector of two colours for the bars
#' @param nBars A numeric fixing the width of the plot to accommodate a certain number of bars
#' @param nLowPoints A logical indicating whether or not plot individual points instead of a box-and-whisker when the number of points is less than 4
#' @details
#' None.
#' @keywords phylloR
#' @return None.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

# Barplot of the subset divided by foi
plotCommunityChanges <- function(cds,soi=NULL,cutoff=0.05,rowLabs=NULL,subtitle=NULL,cols=1:2,nBars=NULL,nLowPoints=TRUE){
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

    par(mar=0.1+c(12,4,4,1))

    counts <- t(cds$counts[soi,])
    ra <- t(100*counts/apply(counts,1,sum))
    ra <- ra[(cds$results[soi,]$padj<cutoff) & !is.na(cds$results[soi,]$padj),,drop=F]
    if(nrow(ra)<1){
        plot.new()
        cat("No significant changes to plot\n")
        return()
    }
    if(is.null(nBars)){
        nBars = 2*nrow(ra)
    }
    rowLabs <- rowLabs[(cds$results[soi,]$padj<cutoff) & !is.na(cds$results[soi,]$padj)]

    ras <- split(as.data.frame(t(ra)),cds$meta[,cds$foi],drop=T)
    ras <- lapply(ras,t)
    if(length(ras)<2){
        cat("All samples have the same factor\n")
        return()
    }
    
    stats <- lapply(ras,function(ra) apply(ra,1,function(x) boxplot(x[x>0],plot=F))) # VIOLIN PLOTS?
    zeros <- lapply(ras,function(ra) apply(ra,1,function(x) sum(x==0)))

    plot(1,type="n",xlim=c(0.5,nBars+0.5),ylim=c(1e-3,1e2),log="y",xlab="",ylab="Relative Abundance",xaxt="none",yaxt="none",main=title)
    abline(h=c(1e-2,1e-1,1e0,1e1,1e2),col="grey")
    axis(1,at=2*(1:nrow(ra))-0.5,labels=rowLabs,las=2)
    axis(2,at=c(1e-3,1e-2,1e-1,1e0,1e1,1e2),labels=c("Undetected","0.01%","0.1%","1%","10%","100%"))

    for(i in 1:(2*nrow(ra))){
        s = stats[[2-(i%%2)]][[(i+(i%%2))/2]]
        z = zeros[[2-(i%%2)]][(i+(i%%2))/2]
        if(nLowPoints & ncol(ras[[2-(i%%2)]])<4){
            points(rep(i,ncol(ras[[2-(i%%2)]])),ras[[2-(i%%2)]][(i+(i%%2))/2,],pch=20,col=cols[2-(i%%2)])
        }else{
            rect(i-0.4,s$stats[2],i+0.4,s$stats[4],col=cols[2-(i%%2)])
            segments(rep(i,2),s$stats[c(1,4)],rep(i,2),s$stats[c(2,5)])
            segments(i-0.4,s$stats[3],i+0.4,s$stats[3],lwd=2)
            points(rep(i,length(s$out)),s$out,pch=20,col=cols[2-(i%%2)])
        }
        points(i,1e-3,cex=2*sqrt(z/nBars),col=cols[2-(i%%2)],pch=20)
        text(i,1.3e-3,z)
    }

    for(i in 1:nrow(ra)){
        text((2*i)-0.5,3e-3,paste("FC: x",formatC(2^cds$results[rownames(ra)[i],]$log2FoldChange,digits=3),"\nP: ",formatC(cds$results[rownames(ra)[i],]$padj,digits=3),sep=""),cex=0.5)
    }
}

#' Function to produce a bar or violin plot of an individual community.
#'
#' @param counts A count table in which the rows are strains or OTUs and the columns are samples. The table should include only a single sample group.
#' @param type Determines the plot style, either 'bar', 'violin', 'swarm', 'barswarm' or 'violinswarm'.
#' @param xlabels  An optional vector of strain names, the default is to use the row names of the count table.
#' @param xcols An optional vector of strain colours, the default is to use rainbow(). 
#' @param res The resolution of the histograms used to create the violins for that plot style.
#' @details
#' For each strain, samples in which they were not observed are counted and separated from the main bar or violin.
#' As a result, the number of samples in each bar or violin varies, and the histogram for each violin is therefore normalised.
#' @keywords phylloR
#' @return A list of the statistics for each strain as generated by boxplot().
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

plotCommunity <- function(counts,type="bar",xlabels=NULL,xcols=NULL,res=50){
    if(!type%in%c("bar","violin","points","swarm","barswarm","violinswarm")){
        cat("Not a valid type of community plot\n")
        return()
    }
    if(is.null(xlabels)){
        xlabels <- rownames(counts)
    }
    if(is.null(xcols)){
        xcols=rainbow(nrow(counts))
    }

    par(mar=0.1+c(12,4,1,1))

    counts <- t(counts)
    rowSums <- apply(counts,1,sum)
    ncts <- 100*counts/rowSums
    medians <- apply(ncts,2,function(x) median(x[x>0]))
    ncts <- ncts[,order(medians,decreasing=T)]
    xlabels = xlabels[order(medians,decreasing=T)]
    xcols = xcols[order(medians,decreasing=T)]

    brks = 10^seq(-2.8,2,length=res)

    stats <- apply(ncts,2,function(x) boxplot(x[x>0],plot=F))
    hists <- apply(ncts,2,function(x) hist(x[x>0],breaks=brks,plot=F))
    zeros <- apply(ncts,2,function(x) sum(x==0))

    plot(1,type="n",xlim=c(1,ncol(ncts)),ylim=c(1e-3,1e2),log="y",xlab="",ylab="Relative Abundance",xaxt="none",yaxt="none")
    abline(h=c(1e-2,1e-1,1e0,1e1,1e2),col="grey")
    axis(1,at=1:ncol(ncts),labels=xlabels,las=2)
    axis(2,at=c(1e-3,1e-2,1e-1,1e0,1e1,1e2),labels=c("Undetected","0.01%","0.1%","1%","10%","100%"))

    for(i in 1:ncol(ncts)){
        if(type%in%c("bar","barswarm")){
            s = stats[[i]]
            if(type=="bar"){
                points(rep(i,length(s$out)),s$out,pch=20,col=xcols[i])
                rect(i-0.4,s$stats[2],i+0.4,s$stats[4],col=xcols[i])
            }else{
                points(rep(i,length(s$out)),s$out,pch=20)
                rect(i-0.4,s$stats[2],i+0.4,s$stats[4])
            }
            segments(rep(i,2),s$stats[c(1,4)],rep(i,2),s$stats[c(2,5)])
            segments(i-0.4,s$stats[3],i+0.4,s$stats[3],lwd=2)
        }else if(type%in%c("violin","violinswarm")){
            h = hists[[i]]
            if(sum(h$counts)>0){
                ycoords = approx(h$mids,n=5*res)$y
                lo <- loess(h$counts/max(h$counts*2)~h$mids,span=0.25)
                pred <- predict(lo,ycoords)
                rpred <- predict(lo,rev(ycoords))
                pred[pred<0] <- 0
                rpred[rpred<0] <- 0
                if(type=="violin"){
                    polygon(c(i,i-pred,i,i+rpred),c(10^-2.8,ycoords,1e2,rev(ycoords)),col=xcols[i])
                }else{
                    polygon(c(i,i-pred,i,i+rpred),c(10^-2.8,ycoords,1e2,rev(ycoords)))
                }
                segments(i-0.4,stats[[i]]$stats[3],i+0.4,stats[[i]]$stats[3],lwd=2)
            }
        }else if(type=="points"){
            points(rep(i,nrow(ncts)),ncts[,i],col=xcols[i],pch=20)
        }else if(type=="swarm"){
            lines(c(i,i),c(2e-3,1e2),col=1,lwd=1)
        }
        points(i,1e-3,cex=4*zeros[i]/nrow(ncts),col=xcols[i],pch=20)
        text(i,1.3e-3,zeros[i])
    }
    if(type%in%c("swarm","barswarm","violinswarm")){
        beeswarm(x=as.data.frame(ncts),add=TRUE,col=xcols,pch=20,corral="wrap")
    }

    return(list(propCounts=ncts,stats=stats))
}

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
plotCommunityPie <- function(counts,strainTaxa,cols,taxLabels=NULL,sort.tax=FALSE,drop=TRUE,...){
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


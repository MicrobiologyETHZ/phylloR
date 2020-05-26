

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

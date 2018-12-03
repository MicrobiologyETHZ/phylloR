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
#' The lists for the include and exclude arguments should consist of vectors named to match the columns of the metadata being considered and contain the valid (or invalid) items.
#' For instance: include=list(initial_condition=c("Control","Test"),time=c("t1","t2")).
#' Samples are filtered first by the include list, then the exclude list. You can of course slice the data set in advance according to more complex criteria.
#' @keywords
#' @return A list containing the original dataset, the DESeq2-analysed dataset and the accompanying results.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

as.ultrametric <- function(phylo){
    # A function to add length to tips so that a tree is ultrametric
    if(!is.rooted(phylo)){
        stop("phylo object is not rooted!")
    }
    rtt <- diag(vcv(phylo))
    extra <- max(rtt)-rtt
    isTip <- phylo$edge[,2] <= Ntip(phylo)
    phylo$edge.length[isTip] <- phylo$edge.length[isTip]+extra
    return(phylo)
}

tipOrder <- function(phylo){
    # A function to find the real order of tips in a phylo object
    isTip <- phylo$edge[,2] <= Ntip(phylo)
    ord <- phylo$edge[isTip,2]
    return(ord)
}

reorderTips <- function(phylo){
    # A function to reorder tip labels THAT DOESN'T WORK
    isTip <- phylo$edge[,2] <= Ntip(phylo)
    ord <- order(phylo$edge[isTip,2])
    phylo$edge[isTip,] <- phylo$edge[isTip,][ord,]
    phylo$edge.length[isTip] <- phylo$edge.length[isTip][ord]
    phylo$tip.label <- phylo$tip.label[ord]
    return(phylo)
}

draw.phylo <- function(x1,y1,x2,y2,phylo,direction="r"){
    # A function to draw a tree of class 'phylo' to fit inside an arbitrary space on an existing plot
    if(direction%in%c("r","u")){
        nodex <- node.depth.edgelength(phylo)
        nodey <- node.height(phylo)
        nodey <- nodey-min(nodey) # because the first branch is always offset for some reason
    }else if(direction%in%c("l","d")){
        nodex <- node.depth.edgelength(phylo)
        nodex <- max(nodex)-nodex
        nodey <- node.height(phylo)
        nodey <- nodey-min(nodey)
    }else{
        stop("Direction not recognised, must be \"r\", \"l\", \"u\" or \"d\"!")
    }

    edgex1 <- c(nodex[phylo$edge[,1]],nodex[phylo$edge[,1]])
    edgey1 <- c(nodey[phylo$edge[,1]],nodey[phylo$edge[,2]])
    edgex2 <- c(nodex[phylo$edge[,1]],nodex[phylo$edge[,2]])
    edgey2 <- c(nodey[phylo$edge[,2]],nodey[phylo$edge[,2]])

    xscale <- (x2-x1)/max(edgex1,edgex2)
    yscale <- (y2-y1)/max(edgey1,edgey2)

    if(direction%in%c("r","l")){
        xscale <- (x2-x1)/max(edgex1,edgex2)
        yscale <- (y2-y1)/max(edgey1,edgey2)

        plotx1 <- x1+(xscale*edgex1)
        ploty1 <- y1+(yscale*edgey1)
        plotx2 <- x1+(xscale*edgex2)
        ploty2 <- y1+(yscale*edgey2)
    }else{
        xscale <- (x2-x1)/max(edgey1,edgey2)
        yscale <- (y2-y1)/max(edgex1,edgex2)

        plotx1 <- x1+(xscale*edgey1)
        ploty1 <- y1+(yscale*edgex1)
        plotx2 <- x1+(xscale*edgey2)
        ploty2 <- y1+(yscale*edgex2)
    }

    segments(plotx1,ploty1,plotx2,ploty2)
}

treatmap <- function(tree,mat,mask=NULL,mask.color="lightgrey",overlay=NULL,aspect.ratio=1/exp(1),tip.labels=NULL,tip.colors=NULL,z.cols=NULL,z.res=NULL,tip.label.width=NULL,mat.label.height=NULL,scale=T,mat.label.color="black",mat.col.order=NULL,mat.hclust=FALSE,...){
    # A function to plot a heatmap beside a phylogenetic tree
    # Check if the matrix is one- or two-sided
    if(sign(min(mat,na.rm=T))==sign(max(mat,na.rm=T))){
        matType <- 1
    }else{
        matType <- 2
    }
    # If no resolution is given, take it from the colors given, or set to default based on matrix sidedness
    if(is.null(z.res)){
        if(is.null(z.cols)){
            z.res <- 33-matType
        }else{
            z.res <- length(z.cols)
        }
    }
    # If no colors are given, use defaults based on matrix sidedness
    if(is.null(z.cols)){
        if(matType==1){
            z.cols <- colorRampPalette(c("white","blue"))(z.res+(z.res%%2))
        }else{
            z.cols <- colorRampPalette(c("blue","white","red"))(z.res+(z.res%%2)-1)
        }
    }
    if(length(z.cols)!=z.res){
        stop("Colors chosen do not cover the given resolution!")
    }
    if(is.null(tip.labels)){
        tip.labels <- tree$tip.label
    }

    # Determine tree size in coordinates and rescale
    treeHeight <- Ntip(tree)
    treeWidth <- vcv(tree)[1,1]
    tree$edge.length <- (tree$edge.length*treeHeight*aspect.ratio)/treeWidth
    treeWidth <- treeHeight*aspect.ratio
    matWidth <- ncol(mat)
    
    # Calculate plot dimensions
    if(is.null(tip.label.width)){
        tip.label.width <- 0.2*treeWidth
    }
    if(is.null(mat.label.height)){
        mat.label.height <- 0.4*treeWidth
        if(mat.hclust){
            mat.label.height = mat.label.height*2
        }
    }
    plotWidth <- treeWidth+matWidth+tip.label.width
    plotHeight <- treeHeight+mat.label.height

    # Use mat.col.order preferentially over mat.hclust to hierarchically cluster the matrix columns
    if(!is.null(mat.col.order)){
        hc = list(order=mat.col.order)
    }else if(mat.hclust){
        d = dist(t(mat))
        hc = hclust(d)
    }else{
        hc = list(order=1:ncol(mat))
    }
    mat <- mat[,hc$order]

    # Determine rearrangement order of matrices to match tree
    treeOrder <- match(tree$tip.label[tipOrder(tree)],rownames(mat))

    # If no mask is given, make one
    if(is.null(mask)){
        mask <- matrix(TRUE,nrow(mat),ncol(mat))
    }else{
        mask <- mask[treeOrder,hc$order]
    }

    # Determine the color scale for the heatmap
    mat <- mat[treeOrder,]
    if(any(mat[mask]<0)){
        matMax <- ceiling(max(abs(mat[mask])))
        bins <- seq(-matMax,matMax,length.out=z.res)
        cells <- apply(mat,2,function(x) cut(x,bins,labels=F))
    }else{
        matMax <- ceiling(max(mat[mask]))
        bins <- seq(0:matMax,length.out=z.res)
        cells <- apply(mat,2,function(x) cut(x,bins,labels=F))
    }

    # Plot the tree
    if(is.null(tip.colors)){
        tip.colors=rainbow(Ntip(tree))
    }
    plot(tree,show.tip.label=F,x.lim=c(0,plotWidth),y.lim=c(-2*as.numeric(scale),plotHeight),...)
    text(treeWidth,pos=4,1:Ntip(tree),tip.labels,col=tip.colors)

    # Plot the heatmap
    mat.colors = z.cols[cells]
    mat.colors[!mask] = mask.color
    rect(rep(treeWidth+tip.label.width+(0:(ncol(mat)-1)),each=nrow(mat)),0.5+rep(0:(nrow(mat)-1),ncol(mat)),rep(treeWidth+tip.label.width+(1:ncol(mat)),each=nrow(mat)),0.5+rep(1:nrow(mat),ncol(mat)),col=mat.colors)
    text(treeWidth+tip.label.width+(0.5+0:(ncol(mat)-1)),nrow(mat)+1,colnames(mat),adj=c(0,0.5),srt=90,col=mat.label.color[hc$order])

    # Add the overlay if requested
    if(!is.null(overlay)){
        overlay <- overlay[treeOrder,hc$order]
        overlay[!mask] <- ""
        text(rep(0.5+treeWidth+tip.label.width+(0:(ncol(mat)-1)),each=nrow(mat)),1+rep(0:(nrow(mat)-1),ncol(mat)),overlay,adj=c(0.5,0.5),cex=0.5)
    }

    # Add a scale bar if requested
    if(scale){
        gradient.rect(treeWidth+tip.label.width,-1,treeWidth+tip.label.width+ncol(mat),0,col=z.cols)
        if(matType==1){
            text(treeWidth+tip.label.width+seq(0,ncol(mat),length.out=9),-1.5,seq(0,matMax,length.out=9),adj=c(0.5,0.5),cex=0.5)
        }else{
            text(treeWidth+tip.label.width+seq(0,ncol(mat),length.out=9),-1.5,seq(-matMax,matMax,length.out=9),adj=c(0.5,0.5),cex=0.75)
        }
    }

    # Add a dendrogram if columns were clustered
    if(mat.hclust){
        draw.phylo(0.5+treeWidth+tip.label.width,nrow(mat)+mat.label.height/2,0.5+treeWidth+tip.label.width+ncol(mat)-1,nrow(mat)+mat.label.height,as.phylo(hc),direction="d")
    }
}


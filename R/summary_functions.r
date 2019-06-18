#' Function to plot a heatmap aligned with a phylogenetic tree.
#'
#' @param phylo An object of class "phylo".
#' @param mat A matrix of data for the heatmap with as many columns as tips in the phylogenetic tree.
#' @param mask An optional matrix of logicals equal in dimension to "mat" indicating whether some cells should not have their values plotted in the heatmap.
#' @param mask.color A color for the masked cells.
#' @param overlay An optional matrix of character strings to display on the cells of the heatmap, for instance significance labels.
#' @param aspect.ratio The aspect ratio width/height in which to plot the phylogenetic tree.
#' @param tip.labels A vector of character strings for tip labels, defaulting to those of the "phylo" object if not provided.
#' @param tip.colors A single color or vector of colors for the tip labels.
#' @param z.cols A vector of colors for the z-axis of the heatmap, see details for the default.
#' @param z.res A resolution for the color vector of the z-axis if automatically calculated.
#' @param tip.label.width A numeric indicating how much space to put between the heatmap and tree for plotting the labels, defaulting to 20\% of the plotted tree width.
#' @param scalebar A logical indicating whether or not to include a scalebar in the plot.
#' @param mat.labels A vector of strings to label the columns of the heatmap; defaults to the column names of "mat".
#' @param mat.label.height A numeric indicating how much space to leave above the heatmap for column labels.
#' @param mat.label.colors A single color or vector of colors for the column labels of the heatmap.
#' @param mat.col.order A numeric vector indicating the order in which to plot the columns of the heatmap; ignored if mat.phylo is provided.
#' @param mat.hclust A logical indicating whether or not to hierarchically cluster the columns of the heatmap; ignored if mat.col.order or mat.phylo is provided.
#' @param mat.phylo An object of class "phylo" which is used to order the columns of the heatmap.
#' @details
#' The default color vector for the z-axis is generated automatically if not provided. Firstly it is determined whether "mat" contains both positive and negative values. If it is one-sided then a vector of 32 colors is made with "colorRampPalette", between white and blue. If it is two-sided, a vector of 31 colors is made between blue, white and red. The argument "z.res" can be used to fix the number of colors in the vector, though it will still generate an odd number if the heatmap is two-sided. Providing "z.cols" will ignore all these determinations and split the color vector evenly between the heatmap's minimum and maximum.
#' @keywords None
#' @return None
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

treatmap <- function(phylo,mat,mask=NULL,mask.color="lightgrey",overlay=NULL,aspect.ratio=1/exp(1),tip.labels=NULL,tip.colors=NULL,z.cols=NULL,z.res=NULL,tip.label.width=NULL,scalebar=TRUE,mat.labels=NULL,mat.label.height=NULL,mat.label.colors="black",mat.col.order=NULL,mat.hclust=FALSE,mat.phylo=NULL,...){
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
    # If no labels are given, default to existing labels
    if(is.null(tip.labels)){
        tip.labels <- phylo$tip.label
    }
    if(is.null(mat.labels)){
        if(!is.null(mat.phylo)){
            mat.labels <- mat.phylo$tip.label
        }else{
            mat.labels <- colnames(mat)
        }
    }

    # Determine tree size in coordinates and rescale
    phyloHeight <- Ntip(phylo)
    phyloWidth <- vcv(phylo)[1,1]
    phylo$edge.length <- (phylo$edge.length*phyloHeight*aspect.ratio)/phyloWidth
    phyloWidth <- phyloHeight*aspect.ratio
    matWidth <- ncol(mat)
    
    # Calculate plot dimensions
    if(is.null(tip.label.width)){
        tip.label.width <- 0.2*phyloWidth
    }
    if(is.null(mat.label.height)){
        mat.label.height <- 0.4*phyloWidth
        if(mat.hclust){
            mat.label.height = mat.label.height*2
        }
    }
    plotWidth <- phyloWidth+matWidth+tip.label.width
    plotHeight <- phyloHeight+mat.label.height

    # Use mat.phylo if given, otherwise mat.col.order preferentially over mat.hclust to hierarchically cluster the matrix columns
    if(!is.null(mat.phylo)){
        hc = list(order=match(mat.phylo$tip.label,colnames(mat)))
    }else if(!is.null(mat.col.order)){
        hc = list(order=mat.col.order)
    }else if(mat.hclust){
        d = dist(t(mat))
        hc = hclust(d)
    }else{
        hc = list(order=1:ncol(mat))
    }
    mat <- mat[,hc$order]

    # Determine rearrangement order of matrices to match tree
    phyloOrder <- match(phylo$tip.label[tipOrder(phylo)],rownames(mat))

    # If no mask is given, make one
    if(is.null(mask)){
        mask <- matrix(TRUE,nrow(mat),ncol(mat))
    }else{
        mask <- mask[phyloOrder,hc$order]
    }

    # Determine the color scale for the heatmap
    mat <- mat[phyloOrder,]
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
        tip.colors=rainbow(Ntip(phylo))
    }
    plot(phylo,show.tip.label=F,x.lim=c(0,plotWidth),y.lim=c(-2*as.numeric(scalebar),plotHeight),...)
    text(phyloWidth,pos=4,1:Ntip(phylo),tip.labels,col=tip.colors)

    # Plot the heatmap
    mat.colors = z.cols[cells]
    mat.colors[!mask] = mask.color
    rect(rep(phyloWidth+tip.label.width+(0:(ncol(mat)-1)),each=nrow(mat)),0.5+rep(0:(nrow(mat)-1),ncol(mat)),rep(phyloWidth+tip.label.width+(1:ncol(mat)),each=nrow(mat)),0.5+rep(1:nrow(mat),ncol(mat)),col=mat.colors)
    text(phyloWidth+tip.label.width+(0.5+0:(ncol(mat)-1)),nrow(mat)+1,mat.labels[hc$order],adj=c(0,0.5),srt=90,col=mat.label.colors[hc$order])

    # Add the overlay if requested
    if(!is.null(overlay)){
        overlay <- overlay[phyloOrder,hc$order]
        overlay[!mask] <- ""
        text(rep(0.5+phyloWidth+tip.label.width+(0:(ncol(mat)-1)),each=nrow(mat)),1+rep(0:(nrow(mat)-1),ncol(mat)),overlay,adj=c(0.5,0.5),cex=0.5)
    }

    # Add a scale bar if requested
    if(scalebar){
        gradient.rect(phyloWidth+tip.label.width,-1,phyloWidth+tip.label.width+ncol(mat),0,col=z.cols)
        if(matType==1){
            text(phyloWidth+tip.label.width+seq(0,ncol(mat),length.out=9),-1.5,seq(0,matMax,length.out=9),adj=c(0.5,0.5),cex=0.5)
        }else{
            text(phyloWidth+tip.label.width+seq(0,ncol(mat),length.out=9),-1.5,seq(-matMax,matMax,length.out=9),adj=c(0.5,0.5),cex=0.75)
        }
    }

    # Add a dendrogram if mat.phylo was provided or columns were clustered
    if(!is.null(mat.phylo)){
        draw.phylo(0.5+phyloWidth+tip.label.width,nrow(mat)+mat.label.height/2,0.5+phyloWidth+tip.label.width+ncol(mat)-1,nrow(mat)+mat.label.height,mat.phylo,direction="d")
    }else if(mat.hclust){
        draw.phylo(0.5+phyloWidth+tip.label.width,nrow(mat)+mat.label.height/2,0.5+phyloWidth+tip.label.width+ncol(mat)-1,nrow(mat)+mat.label.height,as.phylo(hc),direction="d")
    }
}

#' Function to summarise the results of multiple phyllosphere experiments
#' @param cdsList A list object containing compatible data sets created by "makeCDS". The names of each item in the list will be used to label the summary tables. The results in the data sets must contain the same rows in the same order, i.e.: this function is designed to summarise peturbations to a core experiment.
#' @details
#' None
#' @keywords None
#' @return A list containing two data frames: fcMatrix is the summary of log2FoldChange results, pvMatrix is the summary of adjusted p-value results (padj).
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

summariseResults <- function(cdsList){
    resultsList <- lapply(cdsList,function(x) x$results)

    fcMatrix <- lapply(resultsList,function(x) x[,"log2FoldChange",F])
    fcMatrix <- do.call(cbind,fcMatrix)
    colnames(fcMatrix) <- names(resultsList)
    fcMatrix[is.na(fcMatrix)] <- 0

    pvMatrix <- lapply(resultsList,function(x) x[,"padj",F])
    pvMatrix <- do.call(cbind,pvMatrix)
    colnames(pvMatrix) <- names(resultsList)
    pvMatrix[is.na(pvMatrix)] <- 1

    return(list(fcMatrix=fcMatrix,pvMatrix=pvMatrix))
}

#' Function to plot a bipartite graph to summarise strain-strain interactions.
#' @param fcMatrix A matrix of fold changes where the rows are labelled with the affected strains and the columns with the affector strains.
#' @param pvMatrix A matrix of p-values matching "fcMatrix".
#' @param leftPhylo An optional phylo object describing the phylogeny of the affector strains.
#' @param rightPhylo An optional phylo object describing the affected strains.
#' @param leftOrder An optional vector of character strings to place the affector strains in a particular order.
#' @param rightOrder An optional vector of character strings to place the affected strains in a particular order.
#' @param leftLabs An optional vector of character strings to label the affector strains.
#' @param rightLabs An optional vector of character strings to label the affected strains.
#' @param leftCols An optional vector of colors for the affector strain labels.
#' @param rightCols An optional vector of colors for the affected strain labels.
#' @param experiment.type A character string that is either "removal" or "addition" to describe the nature of the perturbations to the core experiment.
#' @param tip.label.width A numeric to indicate how much space should be placed between the tips of the tree (if provided) and graph lines.
#' @param cutoff A numeric to specify a p-value cutoff above which lines will not be plotted.
#' @details
#' The design of the bipartite graph is to show the affector strains on the left, with lines connecting them to affected strains on the right.
#' @keywords None
#' @return None
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

plotBipartiteSummary <- function(fcMatrix,pvMatrix,leftPhylo=NULL,rightPhylo=NULL,leftOrder=NULL,rightOrder=NULL,leftLabs=NULL,rightLabs=NULL,leftCols=NULL,rightCols=NULL,leftScale=1,rightScale=1,experiment.type="removal",tip.label.width=0.1,cutoff=0.05){
    if(is.null(leftOrder)){
        if(!is.null(leftPhylo)){
            leftOrder = leftPhylo$tip.label[leftPhylo$tip.label%in%colnames(fcMatrix)]
        }else{
            leftOrder = colnames(fcMatrix)
        }
    }
    if(is.null(rightOrder)){
        if(!is.null(rightPhylo)){
            rightOrder = rightPhylo$tip.label[rightPhylo$tip.label%in%rownames(fcMatrix)]
        }else{
            rightOrder = rownames(fcMatrix)
        }
    }

    if(is.null(leftLabs)){
        leftLabs <- colnames(fcMatrix)
    }
    if(is.null(rightLabs)){
        rightLabs <- rownames(fcMatrix)
    }
    names(leftLabs) <- colnames(fcMatrix)
    names(rightLabs) <- rownames(fcMatrix)


    if(!experiment.type%in%c("removal","addition")){
        stop("Experiment type must be one of \"removal\" or \"addition\"")
    }
    if(experiment.type=="removal"){
        fcMatrix <- -fcMatrix
    }

    if(is.null(leftCols)){
        leftCols=rep("black",ncol(fcMatrix))
    }
    if(is.null(rightCols)){
        rightCols=rep("black",nrow(fcMatrix))
    }
    names(leftCols) <- colnames(fcMatrix)
    names(rightCols) <- rownames(fcMatrix)

    fcMatrix <- fcMatrix[rightOrder,leftOrder]
    pvMatrix <- pvMatrix[rightOrder,leftOrder]

    height = max(leftScale*ncol(fcMatrix),rightScale*nrow(fcMatrix))

    par(mar=c(0,0,0,0)+0.1)
    plot.new()
    plot.window(xlim=c(-1,1),ylim=c(1,height))

    if(!is.null(leftPhylo)){
        lyoffset <- (height-(leftScale*Ntip(leftPhylo)))/2
        draw.phylo(-1,leftScale+lyoffset,-0.75,(leftScale*Ntip(leftPhylo))+lyoffset,leftPhylo,direction="r",show.tip.label=F)
        text(-0.75,(leftScale*1:ncol(fcMatrix))+lyoffset,leftLabs[leftOrder],pos=4,col=leftCols[leftOrder])
    }else{
        lyoffset = max(0,(leftScale*nrow(fcMatrix))-(rightScale*ncol(fcMatrix)))/2
        text(-0.75,(leftScale*1:ncol(fcMatrix))+lyoffset,leftLabs[leftOrder],pos=2,col=leftCols[leftOrder])
    }
    if(!is.null(rightPhylo)){
        ryoffset <- (height-(rightScale*Ntip(rightPhylo)))/2
        draw.phylo(0.75,rightScale+ryoffset,1,(rightScale*Ntip(rightPhylo))+ryoffset,rightPhylo,direction="l",show.tip.label=F)
        text(0.75,(rightScale*1:nrow(fcMatrix))+ryoffset,rightLabs[rightOrder],pos=2,col=rightCols[rightOrder])
    }else{
        ryoffset = max(0,(leftScale*ncol(fcMatrix))-(rightScale*nrow(fcMatrix)))/2
        text(0.75,(rightScale*1:nrow(fcMatrix))+ryoffset,rightLabs[rightOrder],pos=4,col=rightCols[rightOrder])
    }

    t = seq(0,1,length.out=101)
    for(i in 1:nrow(fcMatrix)){
        for(j in 1:ncol(fcMatrix)){
            s = which(rownames(fcMatrix)==colnames(fcMatrix)[j])
            if(length(s)==0){
                s = 0
            }
            p = matrix(c(-0.75+tip.label.width,0,0,0.75-tip.label.width,(leftScale*j)+lyoffset,(leftScale*j)+lyoffset,(rightScale*i)+ryoffset,(rightScale*i)+ryoffset),ncol=2)
            if((i!=s) & (pvMatrix[i,j]<cutoff)){
                lines(bezier(t,p),col=paste(c("#313695","#FFFFFF","#a50026")[sign(fcMatrix[i,j])+2],"77",sep=""),lwd=abs(fcMatrix[i,j]))
            }
        }
    }
}

#' Function to construct an igraph network from fold-change and p-value matrices
#' @param fcMatrix A matrix of fold changes where the rows are labelled with the strains affected and the columns with the strains either added or removed from the core experiment.
#' @param pvMatrix A matrix of p-values matching "fcMatrix".
#' @param cutoff A numeric indicating the p-value significance cutoff.
#' @param type A string, either "removal" or "addition", indicating the sort of experiment being summarised.
#' @details
#' None
#' @keywords None
#' @return An igraph network object.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

igraphFromSummary <- function(fcMatrix,pvMatrix,cutoff,type="removal"){
    if(type=="removal"){
        multiplier = -1
    }else if(type=="addition"){
        multiplier = 1
    }else{
        stop(type," is not a recognised type, should be \"removal\" or \"addition\".")
    }

    edgeList <- list()
    edgeAttrList <- list()
    for(i in 1:nrow(fcMatrix)){
        for(j in 1:ncol(fcMatrix)){
            if((pvMatrix[i,j]<cutoff) & (rownames(fcMatrix)[i]!=colnames(fcMatrix)[j])){
                edgeList[[length(edgeList)+1]] <- c(colnames(fcMatrix)[j],rownames(fcMatrix)[i])
                edgeAttrList[[length(edgeAttrList)+1]] <- c(multiplier*sign(fcMatrix[i,j]),abs(fcMatrix[i,j]),pvMatrix[i,j])
            }
        }
    }
    edgeList <- do.call(rbind,edgeList)
    edgeAttrList <- do.call(rbind,edgeAttrList)
    colnames(edgeAttrList) <- c("Sign","Weight","Significance")

    network <- graph_from_edgelist(edgeList)
    for(attr in colnames(edgeAttrList)){
        edge_attr(network,attr) <- edgeAttrList[,attr]
    }

    # Add vertices for experiments that had no effect
    for(missing in colnames(fcMatrix)[!colnames(fcMatrix)%in%V(network)$name]){
        network <- add_vertices(network,1,name=missing)
    }

    vertex_attr(network,type) <- vertex_attr(network,"name")%in%colnames(fcMatrix)
    vertex_attr(network,"color") <- leafTaxonomy[vertex_attr(network,"name"),"Color"]
    vertex_attr(network,"name") <- leafTaxonomy[vertex_attr(network,"name"),"Name"]

    return(network)
}


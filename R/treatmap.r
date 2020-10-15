#' Function to plot a heatmap aligned with a phylogenetic tree.
#'
#' @param phylo             An object of class "phylo".
#' @param mat               A matrix of data for the heatmap with as many columns as tips in the phylogenetic tree.
#' @param mask              An optional matrix of logicals equal in dimension to "mat" indicating whether some cells should not have their values plotted in the heatmap.
#' @param mask.color        A color for the masked cells.
#' @param overlay           An optional matrix of character strings to display on the cells of the heatmap, for instance significance labels.
#' @param borders.lwd       An optional matrix of line widths for the borders of each cell
#' @param borders.lty       An optional matrix of line types for the borders of each cell
#' @param aspect.ratio      The aspect ratio width/height in which to plot the phylogenetic tree.
#' @param tip.labels        A vector of character strings for tip labels, defaulting to those of the "phylo" object if not provided.
#' @param tip.colors        A single color or vector of colors for the tip labels.
#' @param z.cols            A vector of colors for the z-axis of the heatmap, see details for the default.
#' @param z.res             A resolution for the color vector of the z-axis if automatically calculated.
#' @param z.sym             A logical indicating whether the z-axis should be forced to be symmetrical.
#' @param tip.label.width   A numeric indicating how much space to put between the heatmap and tree for plotting the labels, defaulting to 20\% of the plotted tree width.
#' @param scalebar          A logical indicating whether or not to include a scalebar in the plot.
#' @param mat.labels        A vector of strings to label the columns of the heatmap; defaults to the column names of "mat".
#' @param mat.label.height  A numeric indicating how much space to leave above the heatmap for column labels.
#' @param mat.label.colors  A single color or vector of colors for the column labels of the heatmap.
#' @param mat.col.order     A numeric vector indicating the order in which to plot the columns of the heatmap; ignored if mat.phylo is provided.
#' @param mat.hclust        A logical indicating whether or not to hierarchically cluster the columns of the heatmap; ignored if mat.col.order or mat.phylo is provided.
#' @param mat.phylo         An object of class "phylo" which is used to order the columns of the heatmap.
#' @details
#' The default color vector for the z-axis is generated automatically if not provided. Firstly it is determined whether "mat" contains both positive and negative values. If it is one-sided then a vector of 32 colors is made with "colorRampPalette", between white and blue. If it is two-sided, a vector of 31 colors is made between blue, white and red. The argument "z.res" can be used to fix the number of colors in the vector, though it will still generate an odd number if the heatmap is two-sided. Providing "z.cols" will ignore all these determinations and split the color vector evenly between the heatmap's minimum and maximum. Setting z.sym to TRUE will make the scale symmetrical about zero, if two-sided.
#' @keywords None
#' @return None
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

treatmap <- function(phylo, mat, mask=NULL, mask.color="lightgrey", overlay=NULL, borders.lwd=NULL, borders.lty=NULL, aspect.ratio=1/exp(1), tip.labels=NULL, tip.colors=NULL, z.cols=NULL, z.res=NULL, z.sym=FALSE, tip.label.width=NULL, scalebar=TRUE, mat.labels=NULL, mat.label.height=NULL, mat.label.colors="black", mat.col.order=NULL, mat.hclust=FALSE, mat.phylo=NULL, ...){
    # Check if the matrix is one- or two-sided
    if(sign(min(mat, na.rm=T))==sign(max(mat, na.rm=T))){
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
            z.cols <- colorRampPalette(c("white", "blue"))(z.res+(z.res%%2))
        }else{
            z.cols <- colorRampPalette(c("blue", "white", "red"))(z.res+(z.res%%2)-1)
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
    phyloWidth <- max(vcv(phylo))
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
        hc = list(order=match(mat.phylo$tip.label, colnames(mat)))
    }else if(!is.null(mat.col.order)){
        hc = list(order=mat.col.order)
    }else if(mat.hclust){
        d = dist(t(mat))
        hc = hclust(d)
    }else{
        hc = list(order=1:ncol(mat))
    }
    mat <- mat[, hc$order]

    # Determine rearrangement order of matrices to match tree
    phyloOrder <- match(phylo$tip.label[tipOrder(phylo)], rownames(mat))

    # If no mask is given, make one
    if(is.null(mask)){
        mask <- matrix(TRUE, nrow(mat), ncol(mat))
    }else{
        mask <- mask[phyloOrder, hc$order]
    }

    # Determine the color scale for the heatmap
    mat <- mat[phyloOrder, ]
    if(any(mat[mask]<0)){
        matMax <- ceiling(max(abs(mat[mask])))
        bins <- seq(-matMax, matMax, length.out=z.res)
        cells <- apply(mat, 2, function(x) cut(x, bins, labels=F))
    }else{
        matMax <- ceiling(max(mat[mask]))
        bins <- seq(0:matMax, length.out=z.res)
        cells <- apply(mat, 2, function(x) cut(x, bins, labels=F))
    }

    # Plot the tree
    if(is.null(tip.colors)){
        tip.colors=rainbow(Ntip(phylo))
    }
    plot(phylo, show.tip.label=F, x.lim=c(0, plotWidth), y.lim=c(-2*as.numeric(scalebar), plotHeight), ...)
    text(phyloWidth, pos=4, 1:Ntip(phylo), tip.labels, col=tip.colors)

    # Plot the heatmap
    if(is.null(borders.lwd)){
        borders.lwd = par("lwd")
    }
    if(is.null(borders.lty)){
        borders.lty = par("lty")
    }
    mat.colors = z.cols[cells]
    mat.colors[!mask] = mask.color
    rect(rep(phyloWidth+tip.label.width+(0:(ncol(mat)-1)), each=nrow(mat)), 0.5+rep(0:(nrow(mat)-1), ncol(mat)), rep(phyloWidth+tip.label.width+(1:ncol(mat)), each=nrow(mat)), 0.5+rep(1:nrow(mat), ncol(mat)), col=mat.colors,border=NA)
    text(phyloWidth+tip.label.width+(0.5+0:(ncol(mat)-1)), nrow(mat)+1, mat.labels[hc$order], adj=c(0, 0.5), srt=90, col=mat.label.colors[hc$order])
    # Additional borders added after for correct overlap
    rect(rep(phyloWidth+tip.label.width+(0:(ncol(mat)-1)), each=nrow(mat)), 0.5+rep(0:(nrow(mat)-1), ncol(mat)), rep(phyloWidth+tip.label.width+(1:ncol(mat)), each=nrow(mat)), 0.5+rep(1:nrow(mat), ncol(mat)), col=NA, border="black",lty=borders.lty, lwd=borders.lwd)

    # Add the overlay if requested
    if(!is.null(overlay)){
        overlay <- overlay[phyloOrder, hc$order]
        overlay[!mask] <- ""
        text(rep(0.5+phyloWidth+tip.label.width+(0:(ncol(mat)-1)), each=nrow(mat)), 1+rep(0:(nrow(mat)-1), ncol(mat)), overlay, adj=c(0.5, 0.5), cex=0.5)
    }

    # Add a scale bar if requested
    if(scalebar){
        gradient.rect(phyloWidth+tip.label.width, -1, phyloWidth+tip.label.width+ncol(mat), 0, col=z.cols)
        if(matType==1){
            text(phyloWidth+tip.label.width+seq(0, ncol(mat), length.out=9), -1.5, seq(0, matMax, length.out=9), adj=c(0.5, 0.5), cex=0.5)
        }else{
            text(phyloWidth+tip.label.width+seq(0, ncol(mat), length.out=9), -1.5, seq(-matMax, matMax, length.out=9), adj=c(0.5, 0.5), cex=0.75)
        }
    }

    # Add a dendrogram if mat.phylo was provided or columns were clustered
    if(!is.null(mat.phylo)){
        draw.phylo(0.5+phyloWidth+tip.label.width, nrow(mat)+mat.label.height/2, 0.5+phyloWidth+tip.label.width+ncol(mat)-1, nrow(mat)+mat.label.height, mat.phylo, direction="d")
    }else if(mat.hclust){
        draw.phylo(0.5+phyloWidth+tip.label.width, nrow(mat)+mat.label.height/2, 0.5+phyloWidth+tip.label.width+ncol(mat)-1, nrow(mat)+mat.label.height, as.phylo(hc), direction="d")
    }
}

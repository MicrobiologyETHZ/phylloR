

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

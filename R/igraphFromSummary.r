

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

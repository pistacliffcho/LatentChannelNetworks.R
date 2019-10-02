#' @title Prep Edgelist
#' @description Prepares edgelist for \code{makeLCN}
#' @param edgeCountList An nx3 matrix of undirected edge pairs of the form 
#' n1, n2, count
#' @noRd
prepEdgeCountList = function(edgeCountList){
  edgeCountList = as.matrix(edgeCountList)
  if(ncol(edgeCountList) != 3){ stop("edgeList should have three columns") }
  storage.mode(edgeCountList) = "integer"
  
  min_n = min(edgeCountList[,1:2])
  if(min_n < 1) stop("Node ids must start with 1")
  max_n = max(edgeCountList[,1:2])
  # Breaking up to and from nodes + counts
  # Note: we are only considering undirected graphs now
  n1 = edgeCountList[,1]
  n2 = edgeCountList[,2]
  cnt = edgeCountList[,3]
  
  flat_ind = (n1-1) + (n2-1) * max_n
  flat_ind_rev = (n1-1) * max_n + (n2-1)
  double_inds = c(flat_ind, flat_ind_rev)
  
  double_cnts = c(cnt, cnt)
  splt_cnts = split(double_cnts, double_inds)
  sum_edges = vapply(splt_cnts, function(x) sum(x), 0)
  flat_index = as.integer(names(splt_cnts) ) 
  n1 = floor(flat_index / max_n) + 1
  n2 = flat_index %% max_n + 1
  
  df = data.frame(j_ind = n2, cnt = sum_edges )
  splt_cnts = split(df, n1 ) 
  empty_mat = matrix(nrow = 0, ncol = 3)
  storage.mode(empty_mat) = "integer"
  ans = rep(list( empty_mat ), max_n)
  for(lname in names(splt_cnts) ){
    this_edge_info = as.matrix( splt_cnts[[lname]] )
    ind = as.integer(lname)
    ans[[ind]] = this_edge_info
  }
  return(ans)
}


#' @title Make Ball-Karrer-Newman (BKN) Model
#' @param edgeCountList An nx2 OR nx3 matrix of edges. See details.
#' @param nDims Number of Latent Communities
#' @param unknownEdges An nx2 matrix of edges for which the count is unknown
#' @details `edgeCountList` should be an nx3 matrix, with columns 1 & 2
#' being node ID's and column 3 being the number of edges. 
#' If an nx2 matrix is supplied instead, it is assumed that all 
#' the edge counts in the list are 1. Note that 
#' the graph is considered to be undirected, so `[1,3,4]` and `[3,1,4]`
#' both imply that there are 4 edges between nodes 1 and 3. 
#' If node pairs appear more than once in edgeCountList, edge counts are summed. 
#' @noRd
makeBKN = function(edgeCountList, nDims = 5, unknownEdges = NULL){
  # If no edge counts (i.e. column 3) provided, assume all counts = 1
  if(ncol(edgeCountList) == 2){
    edgeCountList = cbind(edgeCountList, 1)
  }
  # Check data is (now) nx3
  if(ncol(edgeCountList) != 3)
    stop("edgeCountList should have 3 columns. See ?makeBKN")
  # Getting max node id 
  # Need to allow case where unknownEdges includes a node id greater
  # than those in edgeCountList
  maxKnownNode = max(edgeCountList[,1:2])
  if(!is.null(unknownEdges)){ maxUnknownNode = max(unknownEdges) }
  else{ maxUnknownNode = 0 }
  nNodes = max(c(maxUnknownNode, maxKnownNode))
  # Getting the average edge count. Starting value for unknown edges
  unknownMeanEdges = sum(edgeCountList[,3]) / nNodes^2
  # If no unknownEdges, create a list of vectors of length 0
  if(is.null(unknownEdges)){ 
    unknownList = lapply(rep(0, nNodes), numeric) 
  }
  else{
    # Preparing unknown EdgeList
    unknownList = prepEdgeList(unknownEdges, nNodes) 
    unknownEdges = cbind(unknownEdges, unknownMeanEdges)
    colnames(unknownEdges) = colnames(edgeCountList)
  }
  # BKN algorithm requires both known *and* unknown edges be included
  # in edgeCount list. Here we combine edgeCountList with unknownEdges
  aug_edgeCountList = rbind(edgeCountList, unknownEdges)
  preppedCountList = prepEdgeCountList(aug_edgeCountList)
  theta_init = init_pars(nNodes, nDims)
  if(length(preppedCountList) != length(unknownList)) browser()
  ans = BKN$new(preppedCountList, theta_init, unknownList)
  return(ans)
}

emBKN = function(BKN_mod, max_its = 10000, 
                 tol = 10^-4, 
                 pTol = 10^-5, par = F){
  ans = BKN_mod$em(max_its, tol, pTol, par)
  return(ans)
}
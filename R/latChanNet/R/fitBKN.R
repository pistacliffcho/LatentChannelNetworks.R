#' @title Prep Edgelist
#' @description Prepares edgelist for \code{makeLCN}
#' @param edgeCountList An nx3 matrix of undirected edge pairs of the form 
#' n1, n2, count
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
  ans = split(df, n1 ) 
  ans = lapply(ans, as.matrix)
  
  return(ans)
}


#' @title Make Ball-Karrer-Newman (BKN) Model
#' @param edgeCountList A nx3 matrix of edges. See details.
#' @param nDims Number of Latent Communities
#' @details `edgeCountList` should be an nx3 matrix, with columns 1 & 2
#' being node ID's and column 3 being the number of edges. Note that 
#' the graph is considered to be undirected, so `[1,3,4]` and `[3,1,4]`
#' both imply that there are 4 edges between nodes 1 and 3. 
#' If node pairs appear more than once in edgeCountList, edge counts are summed. 
#' @export
makeBKN = function(edgeCountList, nDims = 5){
  preppedCountList = prepEdgeCountList(edgeCountList)
  nRows = length(preppedCountList)
  theta_init = matrix(
    runif(nRows * nDims, max = 1 / sqrt(nDims)), 
    nrow = nRows)
  
  ans = BKN$new(preppedCountList, theta_init)
  return(ans)
}


#' @export
emBKN = function(BKN_mod, max_its = 10000, tol = 0.0001, par = F, type = 1){
  ans = BKN_mod$em(max_its, tol, par, type)
  return(ans)
}
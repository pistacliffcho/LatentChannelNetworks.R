#' @description  Convert edgeList to flat indices
#' @param edgeList nx2 matrix of edges
#' @param max_n Maximum node ID. 
#' @noRd
ij2flat = function(edgeList, max_n, undirected = T){
  if(missing(max_n)){ stop("max_n must be provided") }
  if(undirected){ edgeList = sort_ij(edgeList) }
  i = edgeList[,1]
  j = edgeList[,2]
  flat = i - 1 + (j - 1) * max_n
  return(flat)
}

#' @description Convert flat indices to edgeList
#' @param flat Flat index, as returned by `ij2flat`. Starts at 0
#' @param max_n Maximum node ID
#' @noRd
flat2ij = function(flat, max_n){
  if(missing(max_n)){
    stop("max_n must be included")
  }
  i = flat %% max_n + 1
  j = floor(flat / max_n) + 1
  ans = cbind(i,j)
  colnames(ans) = c("i", "j")
  return(ans)
}

#' @description Sort edgelist s.t. i >= j
#' @noRd
sort_ij = function(edgeList){
  j_isLess = edgeList[,2] < edgeList[,1]
  
  ans = edgeList
  ans[j_isLess,1] = edgeList[j_isLess, 2]
  ans[j_isLess,2] = edgeList[j_isLess, 1]
  colnames(ans) = c("i", "j")
  return(ans)
}

#' @description Return unique undirected edges from edgeList
#' @param edgeList A nx2 matrix of edges
#' @export
unq_edges = function(edgeList){
  max_n = max(edgeList)
  sorted_el = sort_ij(edgeList)
  flat = ij2flat(sorted_el, max_n, undirected = T)
  unq_flat = unique(flat)
  ans = flat2ij(unq_flat, max_n)
  return(ans)
}

#' Return only unique non-selfloop samples from flat_index
unq_nondiag_flat = function(flat, max_n){
  ij = flat2ij(flat, max_n)
  ij = ij[ ij[,1] != ij[,2],,drop = F ]
  ij_sort = sort_ij(ij)
  flat = ij2flat(ij_sort, max_n)
  ans = unique(flat)
  return(ans)
}

#' Check that there are no missing edges in original edge list
checkMissingList = function(obs_edges, missing_edges, max_n){
  flat_edges = ij2flat(obs_edges, max_n, undirected = T)
  flat_missing = ij2flat(missing_edges, max_n, undirected = T)
  if(any(flat_missing %in% flat_edges))
    stop("Missing edges found in observed edges list!")
}

  
#' @description Sample edges NOT in undirected edgeList 
#' @param edgeList nx2 matrix of edges
#' @param n Number of samples
#' @noRd
sample_nonEdges = function(edgeList, n = 100){
  edgeList = unq_edges(edgeList)
  max_n = max(edgeList)
  sort_edgeList = sort_ij(edgeList)
  flat_index = ij2flat(sort_edgeList, max_n)
  max_flat = max_n^2 - 1
  
  flat_sample = sample(1:max_flat, 4 * n, replace = F)
  flat_sample = unq_nondiag_flat(flat_sample, max_n)
  tries = 0
  flat_unobs = flat_sample[!(flat_sample %in% flat_index)]
  while(tries < 10 & length(flat_unobs) < n){
    tries = tries + 1
    flat_sample = sample(1:max_flat, 4 * n, replace = F)
    flat_sample = unq_nondiag_flat(flat_sample, max_n)
    flat_unobs = flat_sample[!(flat_sample %in% flat_index)]
  }
  if(tries == 10 & length(flat_unobs) < n){
    stop("some sort of error: failed to find missing edges in 10 tries!")
  }
  
  flat_use = flat_unobs[1:n]
  ans = flat2ij(flat_use, max_n)
  return(ans)
}

#' @description Randomly mask edge status
#' @param edgeList nx2 matrix of edges
#' @param nEdges Number of edges to mask
#' @param nNotEdges Number of non-edges that will be masked
#' @description Randomly selects `nEdges` to mask, which 
#' will be removed from the observed edgeList and returned
#' as `masked_edges` and `nNotEdges`, 
#' which don't appear in original edgeList but will returned
#' in `masked_nonEdges`
#' @export
random_splitEdges = function(edgeList, nEdges, nNotEdges){
  edgeList = unq_edges(edgeList)
  nRow = nrow(edgeList)
  edge_indices = sample(1:nRow, nEdges, replace = F)
  edge_pairs = edgeList[edge_indices, ]
  reduced_edges = edgeList[-edge_indices, ]
  nonEdge_pairs = sample_nonEdges(edgeList, nNotEdges)
  
  colnames(reduced_edges) = c("i","j")
  colnames(edge_pairs) = c("i", "j")
  colnames(nonEdge_pairs) = c("i","j")
  
  ans = list(obs_edges = reduced_edges, 
             masked_edges = edge_pairs, 
             masked_nonEdges = nonEdge_pairs)
  return(ans)
}

# Extract mean edges from model
meanEdges = function(mod, edgeList){
 c_edgeList = as.matrix(edgeList) - 1 
 isLCN = is(mod, "Rcpp_LCN")
 if(!isLCN){
   if(!is(mod, "Rcpp_BKN")) 
     stop("mod must be LCN or BKN mod")
 }
 ans = numeric(nrow(edgeList))
 for(i in seq_len(nrow(edgeList))){
   e = c_edgeList[i,]
   if(isLCN){ ans[i] = mod$edgeProb(e[1], e[2]) }
   else{ ans[i] = mod$meanEdges(e[1],e[2]) }
 }
 return(ans)
}
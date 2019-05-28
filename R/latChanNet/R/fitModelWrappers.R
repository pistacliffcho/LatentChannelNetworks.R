split_wEmpty = function(n1, n2){
  n12 = c(n1,n2)
  n21 = c(n2,n1)
  rng = range(n12)
  full_list = list()
  all_names = as.character(1:rng[2])
  full_list[all_names] = -1
  
  splt_list = split(n12, n21)
  splt_names = names(splt_list)
  full_list[splt_names] = splt_list
  
  missing_names = setdiff(all_names, splt_names)
  for(i in seq_along(missing_names)){
    this_name = missing_names[i]
    full_list[[this_name]] = integer()
  }
  return(full_list)
}

getUniqueEdges = function(ind, split_list){
  char_ind = as.character(ind)
  edges = split_list[[char_ind]]
  if(is.null(edges)){
    return(integer(0))
  }
  ans = unique(edges)
  return(edges)
}

prepEdgeList = function(edgeList){
  edgeList = as.matrix(edgeList)
  if(ncol(edgeList) != 2){ stop("edgeList should have only two columns") }
  storage.mode(edgeList) = "integer"
  
  min_n = min(edgeList)
  if(min_n < 1) stop("Node ids must start with 1")
  max_n = max(edgeList)
  # Breaking up to and from nodes 
  # Note: we are only considering undirected graphs now
  n1 = edgeList[,1]
  n2 = edgeList[,2]
  
  # Dropping circles
  drop = n1 == n2
  n1 = n1[!drop]
  n2 = n2[!drop]
  
  # Because undirected, including (n1,n2) and (n2,n1) pairs
  n1_new = c(n1, n2)
  n2_new = c(n2, n1)
  n1 = n1_new
  n2 = n2_new

  splt_edges = split(n2, n1)
  ans = rep(list(integer() ), max_n)
  for(lname in names(splt_edges) ){
    unq_edges = unique(splt_edges[[lname]])
    ind = as.integer(lname)
    ans[[ind]] = unq_edges
  }
  return(ans)
}

#' @title Make Latent Channel Network Model
#' @param edgeList A nx2 matrix of edges
#' @param nDims Number of Latent Channels
#' @export
makeLCN = function(edgeList, nDims = 5){
  preppedEdgeList = prepEdgeList(edgeList)
  nRows = length(preppedEdgeList)
  pmat_init = matrix(
    runif(nRows * nDims, max = 1 / sqrt(nDims)), 
    nrow = nRows)
  
  ans = LCN$new(preppedEdgeList, pmat_init)
  return(ans)
}

#' @title Optimization of Latent Channel Network via EM
#' @param LCN_mod LCN model, output from makeLCN
#' @param iters Maximum iterations
#' @param tol Convergence tolerance
#' @param pTol Tolerance for skipping parameter updates
#' @export
emLCN = function(LCN_mod, iters = 10000, 
                 tol = 10^-4, 
                 pTol = 10^-8){
  ans = LCN_mod$cache_em(iters, tol, pTol)
  return(ans)
}
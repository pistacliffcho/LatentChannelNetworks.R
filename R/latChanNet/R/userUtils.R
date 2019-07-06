#' @title Expected Channel Usage
#' @description Compute expected channel usage between two connected nodes.
#' @param i Index of first node
#' @param j Index of second node
#' @param lcn_mod LCN model
#' @details Computes the expected channel usage between two nodes 
#' \strong{conditional on the two nodes sharing an edge.}
#' @export
computeTheta = function(i, j, lcn_mod){
  # Converting R indices to C indices
  c_i = i - 1; c_j = j - 1
  ans = lcn_mod$computeTheta(c_i,c_j)
  return(ans)
}


#' @title Compute expected connections for node
#' @description Computes the estimated expected connections through each latent channel
#' @param i Index of node
#' @param lcn_mod LCN model
#' @export
computeExpConnects = function(i, lcn_mod){
  # Converting R index to C indices
  c_i = i - 1
  ans = lcn_mod$expectedConnections(c_i)
  names(ans) = paste("Channel", seq_along(ans))
  return(ans)
}

#' Get auc for a model from edges/notEdges list
get_auc = function(mod, edges, notEdges){
  all_edges = rbind(edges, notEdges)
  nEdges = nrow(edges)
  nNotEdges = nrow(notEdges)
  preds = meanEdges(mod, all_edges)
  hasEdge = rep(0, nEdges + nNotEdges)
  hasEdge[seq_along(edges[,1])] = 1
  auc = mltools::auc_roc(preds, hasEdge)
  return(auc)
}

#' Get both in sample and out of sample AUC
get_both_auc = function(mod, 
                        out_edges, 
                        out_notEdges, 
                        in_edges, 
                        in_notEdges){
  out_auc = get_auc(mod, out_edges, out_notEdges)
  in_auc = get_auc(mod, in_edges, in_notEdges)
  ans = data.frame(Out = out_auc, In = in_auc)
  return(ans)
}

#' @description Estimate Out-of-Sample AUC
#' @param edgeList nx2 matrix of edges
#' @param models Character vector of models to use
#' @param nEdgesMasked Number of edges to mask
#' @param nNonEdgesMasked Number of non-edges to mask
#' @export
est_auc = function(edgeList, models = c("LCN", "BKN"),
                   nChan = 10,
                   nEdgesMasked = 400, 
                   nNonEdgesMasked = 400){
  
  colnames(edgeList) = c("i", "j")
  
  split_edges = random_splitEdges(edgeList, 
                                  nEdgesMasked, 
                                  nNonEdgesMasked)
  
  obs_edges = split_edges$obs_edges
  unseen_edges = rbind(split_edges$masked_edges, 
                       split_edges$masked_nonEdges)
  ans = NULL
  
  in_edge_sample = 
    obs_edges[sample(1:nrow(obs_edges), nEdgesMasked),]
  in_nonEdge_sample = 
    sample_nonEdges(rbind(edgeList, 
                          split_edges$masked_nonEdges), 
                    nNonEdgesMasked)
  
  out_edges = split_edges$masked_edges
  out_nonEdges = split_edges$masked_nonEdges
  if("LCN" %in% models){
    lcn_mod = makeLCN(obs_edges, nChan, unseen_edges)
    em_res = emLCN(lcn_mod, 10000, type = "ParEM")
    auc_res = get_both_auc(lcn_mod, 
                      out_edges, 
                      out_nonEdges, 
                      in_edge_sample, 
                      in_nonEdge_sample)
    colnames(auc_res) = paste0("LCN_", colnames(auc_res) )
    ans = auc_res
  }
  if("BKN" %in% models){
    bkn_mod = makeBKN(cbind(obs_edges, 1), 
                      nChan, 
                      unseen_edges)
    em_res = emBKN(bkn_mod, 10000, par = T)
    auc_res = get_both_auc(bkn_mod, 
                           out_edges, 
                           out_nonEdges, 
                           in_edge_sample, 
                           in_nonEdge_sample)
    colnames(auc_res) = paste0("BKN_", colnames(auc_res) )
    if("LCN" %in% models)
      ans = cbind(ans, auc_res)
    else
      ans = auc_res
  }
  return(ans)
}
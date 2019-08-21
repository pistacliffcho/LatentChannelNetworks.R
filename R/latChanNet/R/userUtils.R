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
  preds = mod$predict(all_edges[,1],all_edges[,2])
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
    lcn_mod = makeLatentModel(obs_edges, 
                              nChan, model = "LCN",
                              missingList = unseen_edges)
    em_res = lcn_mod$fit()
    auc_res = get_both_auc(lcn_mod, 
                      out_edges, 
                      out_nonEdges, 
                      in_edge_sample, 
                      in_nonEdge_sample)
    colnames(auc_res) = paste0("LCN_", colnames(auc_res) )
    ans = auc_res
  }
  if("BKN" %in% models){
    bkn_mod = makeLatentModel(cbind(obs_edges, 1), 
                      nChan, model = "BKN",
                      missingList = unseen_edges)
    em_res = bkn_mod$fit()
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





#' @export
simBlockLCN = function(nBlocks = 8, nPerBlock = 32, 
                       nSuper = 8,
                       in_pars = c(2,2), 
                       out_pars = c(1,20), 
                       sparse = 0.75){
  nNodes = nBlocks * nPerBlock
  pmat = matrix(rbeta(nNodes * nBlocks, 
                      out_pars[1], out_pars[2]), 
                nrow = nNodes)
  is_sparse = rbinom(nNodes * nBlocks, 
                     1, prob = sparse)
  pmat[is_sparse == 1] = 0
  
  for(i in 1:nBlocks){
    this_block = 1:nPerBlock + (i-1)*nPerBlock
    pmat[this_block, i] = rbeta(nPerBlock, 
                                in_pars[1], in_pars[2])
  }
  super_pmat = matrix(rbeta(nSuper * nBlocks, 
                            in_pars[1], in_pars[2]), 
                      nrow = nSuper)
  pmat = rbind(pmat, super_pmat)
  blockID = c(rep(1:nBlocks, each = nPerBlock), 
              rep("super", nSuper) )
  edges = simLCN(pmat)
  ans = list(edgeList = edges,
             blockID = blockID,
             chanProbs = pmat)
  return(ans)
}

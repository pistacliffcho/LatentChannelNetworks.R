
# Get auc for a model from edges/notEdges list
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

# Get both in sample and out of sample AUC
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

#' @title Estimate Out-of-Sample AUC
#' @param edgeList nx2 matrix of edges
#' @param models Character vector of models to use
#' @param nChan Number of channels to use
#' @param nEdgesMasked Number of edges to mask
#' @param nNonEdgesMasked Number of non-edges to mask
#' @export
est_auc = function(edgeList, models = c("LCN", "BKN"),
                   nChan = 10,
                   nEdgesMasked = 400, 
                   nNonEdgesMasked = 400){
  em_type = "base"
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
    if(em_type == "base"){
      em_res = lcn_mod$fit()
    }
    else if(em_type == "fast"){
      em_res = emLCN(lcn_mod$cmod, 
                     type = "ParEM", 
                     fast_em = T)
    }
    else if(em_type == "hybrid"){
      em_res1 = emLCN(lcn_mod$cmod, 
                      iters = 100,
                     type = "ParEM", 
                     fast_em = F)
      em_res2 = emLCN(lcn_mod$cmod, 
                     type = "ParEM", 
                     fast_em = T)
    }
    else{ stop("em_type not recognized") }
    
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

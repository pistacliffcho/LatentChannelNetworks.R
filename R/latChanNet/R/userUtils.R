
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


#' Compute sizes of channels
#' 
#' @param mod LatClass model
#' @param type How size is defined. Either 'nodes_using' or 'exp_connects'
#' @description Returns the size of each channel
#' @details The size of each channel can be defined in two ways: 
#' number of nodes that have non-zero attachment to a channel 
#' ('nodes_using') *or* 
#' the expected number of connections through a channel a new node 
#' would have if it connected through that channel with probability 1
#' ('exp_connects')
#' We note that the 'exp_connects' metric is a better description of size, 
#' but 'nodes_using' is more intuitive. 
#' @export 
#' @examples 
#' data(email_data)
#' mod = makeLatentModel(email_data$edgeList, 10, 
#'                       metadata = email_data$meta)
#' mod$fit(fast_em = TRUE)
#' 
#' channel_sizes(mod, "exp_connects") 
channel_sizes = function(mod, type = "nodes_using"){
  node_pars = mod$get_pars(node = TRUE, meta = FALSE)$nodes
  if(type == "nodes_using"){
    binary_usage = node_pars > 0
    ans = colSums(binary_usage)
    names(ans) = colnames(node_pars)
  }
  else if(type == "exp_connects"){
    ans = colSums(node_pars)
  }
  else{
    stop("'type' not recognized. Choices are 'nodes_using' or 'exp_connects'")
  }
  return(ans)
}

#' Estimate Channels Nodes Connect Through
#' @param i Node ids
#' @param j Node ids. If left blank, will select *all* edges with i
#' @param model LatClass model
#' @description Estimates probability two nodes 
#' are connected through a channel, *conditional on them having an edge.*
#' 
#' @examples 
#' data("email_data")
#' mod = makeLatentModel(email_data$edgeList, 10, 
#'                       meta = email_data$meta)
#' mod$fit(fast_em = TRUE)
#' 
#' # Checking channel usage for 
#' # first few edges
#' nodes_1 = email_data$edgeList[1:5, 1]
#' nodes_2 = email_data$edgeList[1:5, 2]
#' chan_connect(nodes_1, nodes_2, mod)
#' 
#' # Checking channel usage for all edges 
#' # for pair of nodes with only a few edges
#' chan_connect(i = c(1000, 1001), model = mod)
#' @export
chan_connect = function(i, j = NULL, model){
  if(is.null(j)){
    expanded_i = NULL
    for(this_i in i){
      these_j = c(model$edgeList_nodesOnly[[this_i]], 
                  model$missingList_nodesOnly[[this_i]])
      these_i = rep(this_i, length(these_j))
      expanded_i = c(expanded_i, these_i)
      j = c(j, these_j)
    }
    i = expanded_i
  }
  if(length(i) == 0){ return(NULL) }
  param_mat = model$get_pars(nodes = TRUE, meta = FALSE)$nodes
  ans = chanConnect(i, j, param_mat, model$model)
  colnames(ans) = paste0("Channel ", seq_len(ncol(ans)))
  rownames(ans) = paste0("Edge ", i, ":", j)
  return(ans)
}

#' Subsets channels that are predictive of metadata
#' 
#' @param model LatClass model
#' @param metanames Vector of names of metadata values to predict
#' @param metavars Vector of column names of metadata to predict
#' @param threshold Minimal parameter value to be considered predictive
#' @param sumFun Summary function: suggest either \code{max} or \code{sum}
#' 
#' @description 
#' Returns both a vector of channels that are predictive of at least one
#' of the metadata values of interest and a parameter matrix of channels by 
#' metadata.
#' 
#' @details 
#' \code{metanames} refers to the individual values we might want to predict, 
#' while \code{metavars} is the column names. 
#' 
#' @examples 
#' data(email_data)
#' mod = makeLatentModel(email_data$edgeList, 20,
#'                       meta = email_data$meta)
#' mod$fit(fast_em = TRUE)
#' 
#' # Returns channels that are predictive 
#' # of dpt == 1 or 2
#' predicts_meta(mod, metanames = c("dpt1", "dpt2") ) 
#' # Returns channels that are predictive 
#' # of *any* dpt
#' predicts_meta(mod, metanames = NULL, metavars = "dpt")
#' @export
predicts_meta = function(model, metanames = NULL, 
                         metavars = NULL, threshold = 0.5, 
                         sumFun = max){
  for(v in metavars){
    metanames = c(metanames, model$metalookup[[v]])
  }
  pars = model$get_pars(nodes = FALSE, meta = TRUE)$meta[metanames,,drop=F]
  sum_vals = apply(pars, 2, sumFun)
  keep = sum_vals > threshold

  ans = list(channels = as.numeric(which(keep)), 
       pars = pars[,keep])
  return(ans)
}
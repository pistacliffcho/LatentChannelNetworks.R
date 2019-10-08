#' @title Predictions from LatClass objects
#' @description Predict edge probabilities and categorical metadata
#' @param mod LatClass model
#' @param i node index
#' @param j Either an node index or metadata colname name
#' @param type Should node pairs or cross of all combinations be predicted 
#' @examples 
#' data(email_data)
#' df = data.frame(dpt = email_data$nodeDpt)
#' # Grouping dpt for brevity
#' df$dpt[df$dpt > 4] = "other"
#' 
#' # Building model and fitting
#' mod = makeLatentModel(email_data$edgeList, 
#'                       nDims = 10, 
#'                       metadata = df)
#' mod$fit(fast_em = T)
#'
#' # Predicting edge pairs
#' predict(mod, i = 1:3, j = 4:2)
#' 
#' # Predicting all combinations of i and j
#' predict(mod, i = 1:3, j = 1:3, type = "cross")
#' 
#' # Predicting metadata 
#' predict(mod, i = 1:3, "dpt")
#' @export
predict.LatClass = function(mod, i, j, type = "pairs"){
  if(type == 'pairs'){
    ans = mod$predict(i,j)
    return(ans)
  }
  else if(type == "cross"){
    ans = mod$cmod$crossEdges(i, j)
    rownames(ans) = paste("Node", i)
    colnames(ans) = paste("Node", j)
    return(ans)
  }
  else{
    stop("type must be 'pair' or 'cross'")
  }
}

LatClass = setRefClass("LatClass", 
                       fields = c("cmod", 
                                  "org_edgeList", 
                                  "org_missingList",
                                  "used_edgeList",
                                  "used_missingList",
                                  "metadata", 
                                  "metanames",
                                  "metalookup",
                                  "modtype", 
                                  "max_node"), 
                       methods = c("fit", 
                                   "predict", 
                                   "softMax_prob",
                                   "plot", 
                                   "llk", 
                                   "get_pars")
)

#' @title Make Latent Structure model
#' @param edgeList An matrix edgelist. Can be nx2 (both) or nx3 (BKN only)
#' @param nDim Number of latent dimensions to use
#' @param model Type of model to fit. Options are "LCN" or "BKN"
#' @param missingList A nx2 matrix edgelist of edges for which the value is unknown
#' @param metadata A data.frame with all factors representing metadata
#' @description Make a latent class model. 
#' Can be used for predicting unknown edge status 
#' and unknown metadata.
#' @details 
#' Fits either a Latent Channels Network (LCN), 
#' or the symmetric low-rank Poisson model of 
#' Ball, Karrer and Newman (BKN). 
#' The model assumes an undirected graph. 
#' 
#' If edges are counts, use the BKN model. 
#' The data format for each row is (i,j, count),
#' with i,j as integer IDs starting at 1. 
#' 
#' If edges are binary, either a BKN or LCN model 
#' may be used, 
#' although an LCN model is somewhat more appropriate. 
#' 
#' LCN model:
#' 
#' Clifford Anderson-Bergman, Phan Nguyen, and Jose Cadena Pico.
#' "Latent Channel Networks", submitted 2019
#' 
#' BKN model: 
#' 
#' Brian Ball, Brian Karrer, and Mark EJ Newman. 
#' "Efficient and principled method for detecting 
#' communities in networks." 
#' Physical Review E 84.3 (2011): 036103.
#' 
#' @examples 
#' data(email_data)
#' # Building model with metadata
#' df = data.frame(dpt = email_data$nodeDpt)
#' model = makeLatentModel(email_data$edgeList, 
#'                         10, 
#'                         metadata = df)
#' # Fitting model
#' model$fit()
#' 
#' # Predicting a two edge probabilities
#' predict(model, )
#' 
#' @export
makeLatentModel = function(edgeList, nDims,
                           model = "LCN",
                           missingList = NULL, 
                           metadata = NULL){
  if(missing(nDims)){ stop("nDims must be specified") }
  # Filling in basic fields
  ans = new("LatClass")
  ans$org_edgeList = edgeList
  ans$org_missingList = missingList 
  ans$metadata = metadata
  ans$metanames = NULL
  ans$modtype = model
  max_node = max(edgeList)
  if(!is.null(missingList)){
    max_node = max(c(max_node, max(missingList)))
  }
  ans$max_node = max_node
  
  # Checking that missing edges are not in original edge list
  if(!is.null(missingList)){ 
    checkMissingList(edgeList[,1:2], missingList, max_node)
  }
  
  # Augmenting graph if metadata provided
  if(!is.null(metadata)){
    if(!is.data.frame(metadata)){
      stop("metadata must be dataframe")
    }
    # If we are using BKN model, 
    # need to add count to augmented data
    count = model == "BKN"
    aug_edges = augWithFactors(edgeList, 
                              metadata, 
                              missingList, 
                              count = count)
    edgeList = aug_edges$edges
    missingList = aug_edges$missingEdges
    ans$used_edgeList = edgeList
    ans$used_missingList = missingList
    ans$metanames = aug_edges$metanames
    ans$metalookup = aug_edges$name_list
  }
  
  if(model == "LCN"){
    ans$cmod = makeLCN(edgeList, nDims, missingList)
  }
  else if(model == "BKN"){
    ans$cmod = makeBKN(edgeList, nDims, missingList)
  }
  else{ stop("model not recognized") }
  return(ans)
}

LatClass$methods(
  llk = function(){
    return(cmod$llk())
  }
)

LatClass$methods(
  set_pars = function(pars){
    if(modtype == "LCN"){ ans = cmod$set_pmat(pars) }
    if(modtype == "BKN"){ ans = cmod$set_theta(pars) }
  }
)

LatClass$methods(
  get_pars = function(){
    # Extracting full matrix
    if(modtype == "LCN"){
      all_pars = cmod$get_pmat()
    }
    if(modtype == "BKN"){
      all_pars = cmod$get_theta()
    }
    # Pulling out + naming node-only parameters
    node_par_inds = seq_len(max_node)
    node_pars = all_pars[node_par_inds,]
    rownames(node_pars) = paste("Node", node_par_inds)
    cnames = paste("Channel", seq_len(ncol(all_pars)))
    colnames(node_pars) = cnames
    ans = list(nodes = node_pars)
    # Adding meta data nodes if they are used
    if(!is.null(metanames)){
      meta_pars = all_pars[-node_par_inds, ]
      rownames(meta_pars) = metanames
      ans[["meta"]] = meta_pars
    }
    return(ans)
  }
)

LatClass$methods(
  rand_start = function(){
    samp_pars = get_pars()
    nChans = ncol(samp_pars)
    nNodes = nrow(samp_pars)
    
    new_vals = matrix(runif(nChans * nNodes, max = sqrt(1 / nChans)), 
                      nrow = nNodes)
    if(any(new_vals > 1) ) browser()
    set_pars(new_vals)
  }
)


LatClass$methods(
  fit = function(iters = 10000,
                 par = T, 
                 pTol = 10^-6, 
                 fast_em = F){
    if(modtype == "LCN"){
      alg_type = "EM"
      if(par){ alg_type = "ParEM" }
      emLCN(cmod, iters, type = alg_type, 
            pTol = pTol, fast_em = fast_em)
    }
    else{
      emBKN(cmod, iters, par = par, 
            pTol = pTol)
    }
  }
)

LatClass$methods(
  plot = function(meta_data, minGrpSize = NULL, 
                  xlab = " ", ylab = " ", 
                  prob_cols =  c("black", "white", "orange", "red"), 
                  greater_col = "purple"){
    heatmapLCN(.self, meta_data, minGrpSize = minGrpSize, 
              xlab = xlab, ylab = ylab, prob_cols = prob_cols, 
              greater_col = greater_col)
  }
)

# Compute max probability class for meta categories
LatClass$methods(
  maxprob = function(i, cat){
    all_probs = meta_probs(i, cat)
    max_ind = numeric(length(i))
    for(i in seq_along(max_ind)){
      max_ind[i] = which.max(all_probs[i,])
    }
    ans = metanames[max_ind]
    return(ans)
  }
)

# Compute softmax probabilities for meta categories
LatClass$methods(
  meta_probs = function(i, cat){
    if(!(cat %in% colnames(metadata))){
      stop("meta type not found. See colnames(mod$metadata) for valid options")
    }
    all_cat = metanames[grep(cat, metanames)]
    ans = matrix(nrow = length(i), 
                         ncol = length(all_cat))
    colnames(ans) = all_cat
    for(ind in seq_along(all_cat) ){
      this_cat = all_cat[ind]
      ans[,this_cat] = predict(i, this_cat)
    }
    rsums = rowSums(ans)
    for(i in seq_along(rsums)){
      if(rsums[i] > 0){
        ans[i,] = ans[i,]/rsums[i]
      }
    }
    return(ans)
  }
)


remove_colname = function(fullnames, colnames){
  nchar_col = nchar(colnames)
  ans = substring(fullnames, nchar_col + 1)
  return(ans)
}

LatClass$methods(
  predict_meta = function(i, meta){
    if(length(meta) != 1) stop("meta must be a scalar")
    all_names = metalookup[[meta]]
    name_inds = match(all_names, metanames)
    js = name_inds + max_node
    ans = cmod$crossEdges(i, js)
    
    # softmax standardize
    ans = ans / rowSums(ans)
    # in case any divide by zeros
    ans[is.na(ans)] = 1 / ncol(ans)
    rownames(ans) = paste("Node", i)
    colnames(ans) = remove_colname(all_names, meta)
    return(ans)
  }
)

LatClass$methods(
  predict = function(i, j){
    if(is.character(j)){
      ans = predict_meta(i,j)
      return(ans)
    }
    if(length(j) == 1 & length(i) > 1){
      j = rep(j, length(i))
    }
    if(length(i) == 1 & length(j) > 1){ 
      i = rep(i, length(j))
    }
    if(length(i) != length(j)){
      stop("length(i) != length(j)")
    }
    if(max(i) > max_node) 
      stop("i outside range of indices")
    if(is.numeric(j)){
      if(max(j) > max_node)
        stop("i outside range of indices")
    }
    ans = meanEdges(cmod, cbind(i,j))
    names(ans) = paste0("Nodes ", i, ":", j)
    return(ans)
  }
)



init_pars = function(nNodes, nDims){
  ans = matrix(runif(nNodes * nDims, max = 1/100), 
               nrow = nNodes)
  row_ind = 1:nNodes
  col_ind = (row_ind %% nDims) + 1
  flat_ind = row_ind + (col_ind - 1) * nNodes
  ans[flat_ind] = runif(nNodes)
  return(ans)
}

expNodeConnectMat = function(mod){
  # Extracting only for properly sized matrix
  ans = mod$get_pars() * 0
  nRows = nrow(ans)
  for(i in 1:nRows){
    ans[i,] = computeExpConnects(i, mod$cmod)
  }
  return(ans)
}

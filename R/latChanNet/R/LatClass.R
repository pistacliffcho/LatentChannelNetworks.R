#' @description Predictions from LatClass objects
#' @param mod LatClass model
#' @param i node index
#' @param j Either an edge index or metadata colname name
#' @export
predict.LatClass = function(mod, i, j, type = "edgeprob"){
  if(type == "edgeprob"){ return(mod$predict(i,j)) }
  if(type == "metaprob"){
    if(length(i) != 1) stop("i must be length 1 if type == 'meta'")
    mNames = grep(j, mod$metanames, value = T)
    ans = mod$predict(i, mNames)
    names(ans) = gsub(j, "", mNames)
    ans = ans / sum(ans)
    return(ans)
  }
  if(type == "metamax"){
    probs = predict(mod, i, j, type = "metaprob")
    max_ind = which.max(probs)
    ans = names(probs)[max_ind]
    if(length(ans) == 0) ans = -1
    return(ans)
  }
  else{
    stop("Unrecognized type. Options are 'edgeprob', 'metaprob' or 'metamax'")
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
                                  "modtype", 
                                  "max_node"), 
                       methods = c("fit", 
                                   "predict", 
                                   "softMax_prob",
                                   "plot", 
                                   "llk", 
                                   "get_pars", 
                                   "mult_fit", 
                                   "resize")
)

LatClass$methods(
  resize = function(new_pars){
    cmod$resize(new_pars)
  }
)

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
    if(modtype == "LCN"){
      ans = cmod$get_pmat()
      return(ans)
    }
    if(modtype == "BKN"){
      ans = cmod$get_theta()
      return(ans)
    }
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
  mult_fit = function(nFits = 5, 
                      nInitIts = 500, 
                      nFinalIts = 10000){
    max_llk = -Inf
    for(i in 1:nFits){
      rand_start()
      res = fit(iters = nInitIts)
      this_llk = llk()
      if(this_llk > max_llk){
        max_llk = this_llk
        max_vals = get_pars()
      }
    }
    set_pars(max_vals)
    fit(iters = nFinalIts)
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
    heatmapLCN(cmod, meta_data, minGrpSize = minGrpSize, 
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


LatClass$methods(
  predict = function(i, j, meta_type = NULL){
    if(!is.null(meta_type)){
      if(meta_type == "prob"){
        ans = meta_probs(i,j)
        return(ans)
      }
      else if(meta_type == "class"){
        ans = maxprob(i,j)
        return(ans)
      }
      else{
        stop("meta_type not recognized. Options are 'prob' or 'class'")
      }
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
    if(is.character(j)){
      name_inds = match(j, metanames)
      if(any(is.na(name_inds))){
        stop("invalid metaname. See $metanames for allowed names")
      }
      j = name_inds + max_node
    }
    ans = meanEdges(cmod, cbind(i,j))
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

#' @export
expNodeConnectMat = function(mod){
  # Extracting only for properly sized matrix
  ans = mod$get_pars() * 0
  nRows = nrow(ans)
  for(i in 1:nRows){
    ans[i,] = computeExpConnects(i, mod$cmod)
  }
  return(ans)
}

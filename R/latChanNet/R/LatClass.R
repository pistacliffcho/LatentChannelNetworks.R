LatClass = setRefClass("LatClass", 
                       fields = c("cmod", 
                                  "edgeList", 
                                  "missingList", 
                                  "metadata", 
                                  "metanames",
                                  "modtype", 
                                  "max_node"), 
                       methods = c("fit", 
                                   "predict")
)

#' @export
makeLatentModel = function(edgeList, nDims,
                           model = "LCN",
                           missingList = NULL, 
                           metadata = NULL){
  if(missing(nDims)){ stop("nDims must be specified") }

  # Filling in basic fields
  ans = new("LatClass")
  ans$edgeList = edgeList
  ans$missingList = missingList 
  ans$metadata = metadata
  ans$metanames = NULL
  ans$modtype = model
  max_node = max(edgeList)
  if(!is.null(missingList)){
    max_node = max(c(max_node, max(missingList)))
  }
  ans$max_node = max_node
  
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
  fit = function(iters = 10000, 
                 par = F){
    if(modtype == "LCN"){
      alg_type = "EM"
      if(par){ alg_type = "ParEM" }
      emLCN(cmod, iters, type = alg_type)
    }
    else{
      emBKN(cmod, iters, par = par)
    }
  }
)

LatClass$methods(
  predict = function(i, j){
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
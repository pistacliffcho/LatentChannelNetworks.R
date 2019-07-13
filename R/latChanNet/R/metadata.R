expandFactors = function(edgeList, meta_data){
  meta_data = as.data.frame(meta_data)
  for(i in 1:ncol(meta_data)){
    meta_data[,i] = as.factor(meta_data[,i])
  }
  colNames = colnames(meta_data)
  ans = NULL
  for(cn in colNames){
    frm_txt = paste0("~", cn, " + 0")
    frm = as.formula(frm_txt)
    mf = model.frame(frm, data = meta_data,
                     na.action = na.pass)
    ans = cbind(ans, model.matrix(frm, mf) )
  }
  return(ans)
}
  
getAugNodes = function(edgeList, 
                       meta_data, 
                       addCount = F, 
                       max_node){ 
  expandedData = expandFactors(edgeList, meta_data)
  metanames = colnames(expandedData)
  aug_nodes = which(expandedData == 1, arr.ind = T)
  aug_nodes[,2] = aug_nodes[,2] + max_node
  miss_aug_nodes = which(is.na(expandedData), arr.ind = T)
  miss_aug_nodes[,2] = miss_aug_nodes[,2] + max_node
  cNames = c("i","j")
  if(addCount){ 
    aug_nodes = cbind(aug_nodes, 1)
    miss_aug_nodes = cbind(miss_aug_nodes, 1)
    cNames = c(cNames, "cnt") 
  }
  
  colnames(aug_nodes) = cNames
  colnames(miss_aug_nodes) = cNames
  
  ans = list(edges = aug_nodes, 
             unknown_edges = miss_aug_nodes, 
             metanames = metanames)
  return(ans)
}

augWithFactors = function(edgeList, meta_data, 
                          missingList, 
                          count = F){
  max_node = max(edgeList)
  if(!is.null(missingList))
    max_node = max(c(max_node, max(missingList) ) )
  augNodes = getAugNodes(edgeList, meta_data, 
                         max_node = max_node, 
                         addCount = count)
  metanames = augNodes$metanames
  colnames(edgeList) = colnames(augNodes$edges)
  augEdges = rbind(edgeList, augNodes$edges)
  missingList = rbind(missingList, augNodes$unknown_edges)
  if(nrow(missingList) == 0)
    missingList = NULL
  ans = list(edges = augEdges, 
             missingEdges = missingList[,1:2], 
             metanames = metanames)
  return(ans)
}
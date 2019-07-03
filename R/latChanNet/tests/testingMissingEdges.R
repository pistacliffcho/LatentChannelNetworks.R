library(latChanNet)

setwd("~/Documents/GitHub/LatentChannelNetworks.R/R/examples/emailNetwork")
edgeList = read.table("email-Eu-core.txt") + 1
dpt = read.table("email-Eu-core-department-labels.txt")[,2] + 1

nChan = 10

unq_edges = function(edgeList){
  for(i in 1:nrow(edgeList)){
    edgeList[i,] = sort(edgeList[i,])
  }
  n = max(edgeList)
  flat_inds = edgeList[,1] - 1 + (edgeList[,2] - 1) * n
  flat_inds = unique(flat_inds)
  n1 = (flat_inds %% n) + 1
  n2 = floor(flat_inds / n) + 1
  ans = cbind(n1, n2)
  return(ans)
}

edgeList = unq_edges(edgeList)

i_is_j = function(x, n){
  i = x %% n + 1
  j = floor( x / n ) + 1
  ans = i == j
  return(ans)
}

sample_nonEdges = function(n, edgeList, directed = F){
  i = edgeList[,1]
  j = edgeList[,2]
  max_n = max(c(i,j))
  if(!directed){
    i_new = c(i,j)
    j = c(j,i)
    i = i_new
  }
  flat_index = i - 1 + (j-1) * max_n
  max_flat = max_n^2 - 1
  
  flat_sample = sample(1:max_flat, 4 * n)
  flat_sample = flat_sample[!i_is_j(flat_sample, max_n)]
  tries = 0
  flat_unobs = flat_sample[!(flat_sample %in% flat_index)]
  while(tries < 10 & length(flat_unobs) < n){
    tries = tries + 1
    flat_sample = sample(1:max_flat, 4 * n)
    flat_sample = flat_sample[!i_is_j(flat_sample, max_n)]
    flat_unobs = flat_sample[!(flat_sample %in% flat_index)]
  }
  if(tries == 10 & length(flat_unobs) < n){
    stop("some sort of error: failed to find missing edges in 10 tries!")
  }
  
  flat_use = flat_unobs[1:n]
  new_i = flat_use %% max_n + 1
  new_j = floor( flat_use / max_n ) + 1
  ans = cbind(i = new_i, j = new_j)
  return(ans)
}

random_splitEdges = function(orgList, nEdges, nNotEdges){
  nRow = nrow(orgList)
  edge_indices = sample(1:nRow, nEdges, replace = F)
  
}
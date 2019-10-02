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


getDeg = function(mod){
  theta = mod$get_theta()
  nRow = nrow(theta)
  ans = NULL
  for(i in 1:nRow){
    ans[i] = mod$expectedDegree(i)
  }
  return(ans)
}

edgeList = unq_edges(edgeList)
obsDeg = NULL
for(i in 1:max(edgeList)){
  obsDeg[i] = sum(edgeList == i)
}


scale_model = function(mod, obs_deg){
  tot_degs = sum(obs_deg)
  fit_degs = sum(getDeg(mod))
  rescale_factor = sqrt( tot_degs / fit_degs ) 
  
  old_theta = mod$get_theta()
  new_theta = old_theta * rescale_factor
  mod$set_theta(new_theta)
}

countList = cbind(edgeList, 1)

iter = 100

bkn_mod = makeBKN(countList, nChan)
scale_model(bkn_mod, obsDeg)
emBKN(bkn_mod, max_its = iter, type = 2)
fit_degs = getDeg(bkn_mod)
cat("Fitted degs before rescale =", sum(fit_degs),
    "LLK =", round(bkn_mod$llk(), 3), "\n")
fit_degs = getDeg(bkn_mod)
cat("Fitted degs after rescale =", sum(fit_degs),
    "LLK =", round(bkn_mod$llk(), 3), "\n")

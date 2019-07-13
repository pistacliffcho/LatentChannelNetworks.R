library(latChanNet)
library(RcppParallel)
setThreadOptions(numThreads = 3)

setwd("~/Documents/GitHub/LatentChannelNetworks.R/R/examples/emailNetwork")
edgeList = read.table("email-Eu-core.txt") + 1

nChan = 10

mod = makeLCN(edgeList, nChan)
system.time( res <- emLCN(mod, 10000, type = "ECM") )
res
mod$llk()

mod = makeLCN(edgeList, nChan)
system.time( res <- emLCN(mod, 10000, type = "EM") )
res
mod$llk()


lcn_mod = makeLCN(edgeList, nChan)
system.time( res <- emLCN(lcn_mod, 10000, type = "ParEM") )
res
lcn_mod$llk()


dpt = read.table("email-Eu-core-department-labels.txt")[,2] + 1
#heatmapLCN(mod, dpt, minGrpSize = 15)



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
countList = cbind(edgeList, 1)
set.seed(1)
bkn_mod = makeBKN(countList, nChan)
system.time( res <- emBKN(bkn_mod, par = F, type = 1) )
set.seed(1)
bkn_mod2 = makeBKN(countList, nChan)
system.time( res2 <- emBKN(bkn_mod2, par = F, type = 2) )
bkn_mod2$llk()
pre_t = bkn_mod2$get_theta()
node_llk = NULL
for(i in 0:1004){
  node_llk[i+1] = bkn_mod2$node_llk(i)
}

system.time( res <- emBKN(bkn_mod2, par = F, type = 1) )
node_llk2 = NULL
for(i in 0:1004){
  node_llk2[i+1] = bkn_mod2$node_llk(i)
}

node_llk2 = NULL
heatmapLCN(bkn_mod2, dpt, minGrpSize = 15)

bkn_degs = NULL
for(i in 1:max(edgeList)){
  bkn_degs[i] = bkn_mod2$expectedDegree(i)
}

lcn_degs = NULL
for(i in 1:max(edgeList)){
  lcn_degs[i] = lcn_mod$expectedDegree(i)
}

true_degs = NULL
for(i in 1:max(edgeList)){
  true_degs[i] = sum(edgeList == i)
}

plot(lcn_degs, true_degs)
lines(c(0, 500), c(0, 500), col = 'red')

plot(bkn_degs, true_degs)

heatmapLCN(lcn_mod, dpt)
heatmapLCN(bkn_mod, dpt)


meta_data = data.frame(dpt = dpt)
md = makeLatentModel(edgeList, 10, metadata = meta_data)
md$fit()
md$predict(1, 2)
md$predict(1, "dpt1")
md$predict(2, paste0("dpt", 1:42) )





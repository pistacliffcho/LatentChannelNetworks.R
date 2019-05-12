library(latChanNet)

edgeList = rbind(1:2, 2:3, c(1,3), c(1,4))
pEdgeList = latChanNet:::prepEdgeList(edgeList)
pEdgeList
lcn = makeLCN(edgeList, 2)
lcn$llk()
res1 = lcn$cache_em(10, .01, .01)
lcn$llk()

res2 = emLCN(lcn)
lcn$llk()
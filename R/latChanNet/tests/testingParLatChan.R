library(latChanNet)
library(igraph)
library(RcppParallel)
setThreadOptions(2)

pm <- cbind( c(.1, .01), c(.01, .5) )
g <- sample_sbm(50000, pref.matrix=pm, block.sizes=c(30,70) * 500)
grp = c(rep(1, 15000), rep(2,35000))
el = as_edgelist(g)


mod = makeLCN(el)
mod$llk()
system.time( emLCN(mod, 10) )
mod$llk()
system.time( emLCN(mod, 10, use_par = T) )
mod$llk()

heatmapLCN(mod, grp)

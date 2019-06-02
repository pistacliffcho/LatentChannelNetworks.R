library(latChanNet)
library(igraph)
library(RcppParallel)
setThreadOptions(2)

pm <- cbind( c(.1, .01), c(.01, .5) )
g <- sample_sbm(5000, pref.matrix=pm, block.sizes=c(30,70) * 50)
el = as_edgelist(g)


mod = makeLCN(el)
mod$llk()
system.time( emLCN(mod, 10) )
mod$llk()
system.time( emLCN(mod, 10, use_par = T) )
mod$llk()

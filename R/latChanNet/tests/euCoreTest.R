library(latChanNet)
library(RcppParallel)
setThreadOptions(numThreads = 4)

setwd("~/Documents/GitHub/LatentChannelNetworks/R/examples/emailNetwork")
edgeList = read.table("email-Eu-core.txt") + 1

mod = makeLCN(edgeList, 10)
system.time( res <- emLCN(mod, 10000, type = "ECM") )
res
mod$llk()

mod = makeLCN(edgeList, 10)
system.time( res <- emLCN(mod, 10000, type = "EM") )
res
mod$llk()


mod = makeLCN(edgeList, 10)
system.time( res <- emLCN(mod, 10000, type = "ParEM") )
res
mod$llk()


dpt = read.table("email-Eu-core-department-labels.txt")[,2] + 1
heatmapLCN(mod, dpt, minGrpSize = 15)

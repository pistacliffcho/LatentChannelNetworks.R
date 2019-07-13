library(latChanNet)
library(RcppParallel)
setThreadOptions(numThreads = 3)

setwd("~/Documents/GitHub/LatentChannelNetworks.R/R/examples/emailNetwork")
edgeList = read.table("email-Eu-core.txt") + 1
dpt = read.table("email-Eu-core-department-labels.txt")[,2] + 1

mask = sample(1:1005, 100, replace = F)
dpt_mask = dpt
dpt_mask[mask] = NA


meta_data = data.frame(dpt = as.factor(dpt_mask))
aug_res = latChanNet:::augWithFactors(edgeList, meta_data, NULL)
md = makeLatentModel(edgeList, 40, metadata = meta_data)
md$fit()

fitted_dpt = NULL
for(i in 1:1005){
  fitted_dpt[i] = predict(md, i, "dpt", type = "metamax")
}

is_correct = dpt == fitted_dpt
acc = NULL
acc["inSample"] = mean(is_correct[-mask])
acc["outOfSample"] = mean(is_correct[mask])
acc

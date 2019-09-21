library(latChanNet)

data("email_data")
mod = makeLatentModel(email_data$edgeList, 16)
mod$fit()
expMat = latChanNet:::expNodeConnectMat(mod)
pmat = mod$get_pars()

# These two should be highly correlated
cor(as.numeric(pmat), as.numeric(expMat))
plot(as.numeric(pmat), as.numeric(expMat))

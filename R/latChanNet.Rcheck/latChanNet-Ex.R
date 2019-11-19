pkgname <- "latChanNet"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('latChanNet')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("chan_connect")
### * chan_connect

flush(stderr()); flush(stdout())

### Name: chan_connect
### Title: Estimate Channels Nodes Connect Through
### Aliases: chan_connect

### ** Examples

data("email_data")
mod = makeLatentModel(email_data$edgeList, 10, 
                      meta = email_data$meta)
mod$fit(fast_em = TRUE)

# Checking channel usage for 
# first few edges
nodes_1 = email_data$edgeList[1:5, 1]
nodes_2 = email_data$edgeList[1:5, 2]
chan_connect(nodes_1, nodes_2, mod)

# Checking channel usage for all edges 
# for first two nodes
chan_connect(i = c(1000, 1001), model = mod)



cleanEx()
nameEx("channel_sizes")
### * channel_sizes

flush(stderr()); flush(stdout())

### Name: channel_sizes
### Title: Compute sizes of channels
### Aliases: channel_sizes

### ** Examples

data(email_data)
mod = makeLatentModel(email_data$edgeList, 10, 
                      metadata = email_data$meta)
mod$fit(fast_em = TRUE)

channel_sizes(mod, "exp_connects") 



cleanEx()
nameEx("latChanNet-package")
### * latChanNet-package

flush(stderr()); flush(stdout())

### Name: latChanNet-package
### Title: Latent Network Models for edge and metadata prediction.
### Aliases: latChanNet-package latChanNet
### Keywords: network

### ** Examples

  ## Not run: 
##D      ## Optional simple examples of the most important functions
##D      ## These can be in \dontrun{} and \donttest{} blocks.   
##D   
## End(Not run)



cleanEx()
nameEx("makeLatentModel")
### * makeLatentModel

flush(stderr()); flush(stdout())

### Name: makeLatentModel
### Title: Make Latent Structure model
### Aliases: makeLatentModel

### ** Examples

data(email_data)
# Building model with metadata
model = makeLatentModel(email_data$edgeList, 
                        10, 
                        metadata = email_data$meta)
# Fitting model
model$fit()

# Predicting two edge probabilities
predict(model, i = c(2,3), j = c(4,5))




cleanEx()
nameEx("predict.LatClass")
### * predict.LatClass

flush(stderr()); flush(stdout())

### Name: predict.LatClass
### Title: Predictions from LatClass objects
### Aliases: predict.LatClass

### ** Examples

data(email_data)

# Building model and fitting
mod = makeLatentModel(email_data$edgeList, 
                      nChans = 10, 
                      metadata = email_data$meta)
mod$fit(fast_em = TRUE)

# Predicting edge pairs
predict(mod, i = 1:3, j = 4:2)

# Predicting all combinations of i and j
predict(mod, i = 1:3, j = 1:3, type = "cross")

# Predicting metadata 
# Subsetting for brevity
predict(mod, i = 1:3, "dpt")[,1:5]



cleanEx()
nameEx("predicts_meta")
### * predicts_meta

flush(stderr()); flush(stdout())

### Name: predicts_meta
### Title: Subsets channels that are predictive of metadata
### Aliases: predicts_meta

### ** Examples

data(email_data)
mod = makeLatentModel(email_data$edgeList, 20,
                      meta = email_data$meta)
mod$fit(fast_em = TRUE)

# Returns channels that are predictive 
# of dpt == 1 or 2
predicts_meta(mod, metanames = c("dpt1", "dpt2") ) 
# Returns channels that are predictive 
# of *any* dpt
predicts_meta(mod, metanames = NULL, metavars = "dpt")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

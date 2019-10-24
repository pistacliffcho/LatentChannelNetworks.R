pkgname <- "latChanNet"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('latChanNet')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
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

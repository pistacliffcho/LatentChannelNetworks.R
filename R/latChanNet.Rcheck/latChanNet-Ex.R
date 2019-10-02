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
### Title: A short title line describing what the package does
### Aliases: latChanNet-package latChanNet
### Keywords: package

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
df = data.frame(dpt = email_data$nodeDpt)
model = makeLatentModel(email_data$edgeList, 
                        10, 
                        metadata = df)
# Fitting model
model$fit()

# Predicting a two edge probabilities
predict(model, )




cleanEx()
nameEx("predict.LatClass")
### * predict.LatClass

flush(stderr()); flush(stdout())

### Name: predict.LatClass
### Title: Predictions from LatClass objects
### Aliases: predict.LatClass

### ** Examples

data(email_data)
df = data.frame(dpt = email_data$nodeDpt)
# Grouping dpt for brevity
df$dpt[df$dpt > 5] = "other"
# Building model and fitting
mod = makeLatentModel(email_data$edgeList, 
                      nDims = 10, 
                      metadata = df)
mod$fit(fast_em = T)

# Predicting edge pairs
predict(mod, i = 1:3, j = 1:3)

# Predicting all combinations of i and j
predict(mod, i = 1:3, j = 1:3, type = "cross")

# Predicting meta data 
predict(mod, i = 1:3, "dpt")



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

library(latChanNet)
library(igraph)
library(grid)

set.seed(123)

plot_width = 10
plot_height = 7

## NOTE: TO GET PLOTS TO SAVE CORRECTLY, MUST USED R, NOT RSTUDIO!
## ALSO NEEDED TO RUN TWICE FOR SOME REASON

setwd("/Users/andersonberg1/Documents/GitHub/LatentChannelNetworks/paper")

## Sampling SBM
p_in = 0.25
p_out = 0.025
nGrp = 10
nPerGrp = 100
block_membership = rep(1:nGrp, each = nPerGrp)
# Constructing probability matrix
p_mat = matrix(p_out, nrow = nGrp, ncol = nGrp)
diag(p_mat) = p_in
# Sampling
sbm1 = sample_sbm(nGrp * nPerGrp, p_mat, rep(nPerGrp, nGrp) )
sbm1_edgeList = as_edgelist(sbm1)

# Fitting
mod1 = makeLCN(sbm1_edgeList, nGrp)
emLCN(mod1, iters = 10000)
# Plotting
pdf("baseSBM.pdf", width = plot_width, height = plot_height)
p = heatmapLCN(mod1, grp = block_membership, 
           xlab = "Nodes by Block", 
           ylab = "Channels", 
           main = "")
dev.off()

# Constructing probability matrix
p_mat = matrix(p_out, nrow = nGrp, ncol = nGrp)
diag(p_mat) = p_in


## Sampling SBM with highly connected block
p_mat_aug = cbind(p_mat, p_in)
p_mat_aug = rbind(p_mat_aug, p_in)
block_membership = c(block_membership, 
                     rep(nGrp + 1, nPerGrp))
sbm2 = sample_sbm( (nGrp + 1) * nPerGrp, 
                   p_mat_aug, 
                   rep(nPerGrp, nGrp + 1) )
sbm2_edgeList = as_edgelist(sbm2)


# Fitting
mod2 = makeLCN(sbm2_edgeList, nGrp)
emLCN(mod2, iters = 10000)
# Plotting
heatmapLCN(mod2, grp = block_membership, 
           xlab = "Nodes by Block", 
           ylab = "Channels", 
           main = 
             "Channel Frequency by Node for Stochastic Block Model")

pdf("augSBM.pdf", width = plot_width, height = plot_height)
p = heatmapLCN(mod2, grp = block_membership, 
           xlab = "Nodes by Block", 
           ylab = "Channels", 
           main = 
             "")
dev.off()

# LatentChannelNetworks

R and Julia implementation of the [Latent Channel Networks](https://arxiv.org/abs/1906.04563) EM algorithm. 

The Latent Channel Network (LCN) is a model for undirected graphs in which nodes share an observed edge if they make a connection through *at least* one latent channel. This is meant to reflect modern social dynamics; we are likely to make social connections if we communicate through some group, whether it be due to work, hobbies, familial ties, etc., but it is not necessary to share *all* social groups. 

We can describe this model mathematically by giving each node in the graph a vector of probabilities of length *K*. The probability of two nodes making a connection through the k'th channel is the product of the k'th element of their probability vectors. For here, the probability two nodes share an edge is one minus the probability they make **no** connections through any of the *K* channels. 

In emailExample.md, we provide sample R code for analyzing a classic email network at an European university. Analysis of these results can be found in the paper link above, along with an analysis of a large Facebook network not included in this repository. 

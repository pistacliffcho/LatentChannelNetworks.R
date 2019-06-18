# LatentChannelNetworks

R and Julia implementation of the [Latent Channel Networks](https://arxiv.org/abs/1906.04563) EM algorithm. 

The Latent Channel Network (LCN) is a model for undirected graphs in which nodes share an observed edge if they make a connection through *at least* one latent channel. This is meant to reflect modern social dynamics; we are likely to make social connections if we communicate through some group, whether it be due to work, hobbies, familial ties, etc., but it is not necessary to share *all* social groups. 

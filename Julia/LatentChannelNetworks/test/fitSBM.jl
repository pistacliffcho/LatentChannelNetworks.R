using Revise
using LatentChannelNetworks
using Test

nGrps = 10
nPerGrp = 20
pin = 0.25
pout = 0.05

sbm_edgeList = simSBM(nGrps, nPerGrp, pin, pout)
lcn = makeLatChan(sbm_edgeList, nPerGrp)
starting_llk = computeLLK(lcn)
res = em_cached!(lcn)
finishing_llk = computeLLK(lcn)

@test finishing_llk > starting_llk

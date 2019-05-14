#include("../../src/LatentHubNetworks.jl")
#using Main.LatentHubNetworks
using LatentHubNetworks
using Random
using Plots
gr()
Random.seed!(123)

nGrps = 10
nPerGrp = 100
pIn = 0.25
pOut = 0.025
sbm = simSBM(nGrps, nPerGrp, pIn, pOut)

grpID = []
for i in 1:(nGrps * nPerGrp)
  push!(grpID, 1 + mod(i, nGrps))
end

lmod1 = makeLatHub(sbm, nGrps)
initialize_cache(lmod1, true)

res = em_cached!(lmod1, 5000)
sbm1Plot = heatMap(lmod1, grpID, :magma)
title!("Channel Usage: SBM")
xlabel!("Node")
ylabel!("Channel")
png(sbm1Plot, "baseSBM.png")

sbm2 = augSBM(nPerGrp, pIn, pIn, sbm)
grpID2 = copy(grpID)
for i in 1:nPerGrp
  push!(grpID, nGrps + 1)
end

lmod2 = makeLatHub(sbm2, nGrps)
initialize_cache(lmod2, true)
em_cached!(lmod2, 5000)
sbm2Plot = heatMap(lmod2,grpID, :magma)
title!("Channel Usage: SBM + High Degree Group")
xlabel!("Node")
ylabel!("Channel")
savefig(sbm2Plot, "augSBM.svg")

include("../EMHub.jl")

using Plots
gr()

# Setting seed for replication
Random.seed!(1)

# Control parameters for Stochastic Block Model
nGrps = 5
nPerGrp = 40
pIn = 0.5
pOut = 0.1

smbEdgeList = simSBM(nGrps, nPerGrp, pIn, pOut)

eHub = makeLatHub(smbEdgeList, nGrps)
true_group = repeat( (1:nGrps), nPerGrp)

println("Starting LLK = ", computeLLK(eHub) )
for i in 1:1000
    em_all!(eHub)
end
println("Ending LLK = ", computeLLK(eHub) )

p1 = plotEachDim(eHub, true_group, (1:nGrps) )
png(p1, "sbmExample.png")

nOutliers = 5
sbm_withOutliers = augSBM(nOutliers, 0.5, 0.5, smbEdgeList)

eHub2 = makeLatHub(sbm_withOutliers, nGrps)
println("Starting LLK = ", computeLLK(eHub2) )
for i in 1:1000
    em_all!(eHub2)
end
println("Ending LLK = ", computeLLK(eHub2) )

true_group_aug = cat(true_group,
                     repeat([nGrps + 1], nOutliers),
                     dims = 1 )

p1 = plotEachDim(eHub2, true_group_aug, (1:nGrps) )
png(p1, "sbmWOutliersExample.png")

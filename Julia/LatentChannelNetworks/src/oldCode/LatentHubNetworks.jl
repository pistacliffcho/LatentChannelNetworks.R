module LatentHubNetworks

#function hello_world():
#    println("hello stupid world")
#end
#export hello_world

export EMHub, makeLatHub, probEdge, computeLLK, em_all!, getConditionalHubProb,
    dropLoops, simSBM, augSBM, prepEdges,
    em_cached!, makeTransposeMap, initialize_cache,
    findTransposeInds, getTransposeInds,
    plotEachDim, partionWeight, heatMap,
    organizeWeights, prepWeights4HeatMap,
    allPartionWeights,
    hubSize, computeTheta, hubSizes,
    relabelBySize, expHubConnections


    include("edgeListUtils.jl")
    include("EMHub.jl")
    include("PrecomputedSpeedup.jl")
    include("UserUtils.jl")


end # module

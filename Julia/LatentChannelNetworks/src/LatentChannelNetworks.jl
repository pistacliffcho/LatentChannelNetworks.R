module LatentChannelNetworks

include("edgeListUtils.jl")
include("baseModel.jl")
include("probFunctions.jl")
include("cacheEM.jl")
include("UserUtils.jl")

export simSBM, augSBM
export makeLatChan, probEdge, hubSizes, computeTheta
export em_cached!, computeLLK

end # module

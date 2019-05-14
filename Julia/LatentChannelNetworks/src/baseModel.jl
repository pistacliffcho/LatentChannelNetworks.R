# Latent Channel Class
mutable struct LatentChannelNetwork
    # Number of Nodes in Graph
    nNodes::Int64
    # Number of channels
    dim::Int64
    # List of pairs
    edgeList::Array{Array{Int64, 1}}
    # Matrix of channel frequencies
    pmat::Matrix{Float64}
    # Column mean of pmat
    pbar::Vector{Float64}

    # Cached table of edge probabilities for observed edges
    cache_probs::Vector{Vector{Float64}}
    # Reverse lookup tool so we can easily update i,j and j,i
    cache_map::Vector{Vector{Int64}}
end




## BASIC MODEL QUERIES

# Extracts a boolean vector of whether node is connected to others
function getEdgeBools(i::Int64, lcn::LatentChannelNetwork)::Vector{Bool}
    # Getting number of nodes
    nNodes = lcn.nNodes
    # Creating vector of booleans
    ans = fill(false, nNodes)
    # Extracting edges for this node
    these_edges = lcn.edgeList[i]
    # Setting nodes equal to true
    for j in 1:length(these_edges)
        ans[these_edges[j]] = true
    end
    return(ans)
end

# Compute Log-Likelihood contribution from single node
function computeNodeLLK(i::Int64, lcn::LatentChannelNetwork)::Float64
    nNodes = lcn.nNodes
    ans = 0.0
    these_edges = lcn.edgeList[i]
    hasEdgeVec = getEdgeBools(i, lcn)
    for j in 1:nNodes
        if i == j
            continue
        end
        this_edgeProb = probEdge(i, j, lcn)
        if hasEdgeVec[j]
            ans = ans + log(this_edgeProb)
        else
            ans = ans + log(1.0 - this_edgeProb)
        end
    end
    return(ans)
end


# Computes Log Likelihood
function computeLLK(lcn::LatentChannelNetwork)::Float64
    nNodes = lcn.nNodes
    ans = 0.0
    for i in 1:nNodes
        ans = ans + computeNodeLLK(i, lcn)
    end
    ans = ans / 2.0
    return(ans)
end


# Getter
function get_cache_probs(i::Int64, lcn::LatentChannelNetwork)::Vector{Float64}
    return( lcn.cache_probs[i] )
end
# Setter
function set_cache_probs!(probs::Vector{Float64}, i::Int64, lcn::LatentChannelNetwork)
    updateCacheProb!(probs, i, lcn)
end




# Sparse update tools
function makeTransposeMap(edgeList::Vector{Vector{Int64}})::Vector{Vector{Int64}}
    nNodes = length(edgeList)
    t_map = Vector{Vector{Int64}}(undef, nNodes)
    for i in 1:nNodes
        these_edges = edgeList[i]
        nEdges = length(these_edges)
        this_map = Vector{Int64}(undef, nEdges)
        for j in 1:nEdges
            j_ind = these_edges[j]
            this_map[j] = findTransposeInds(i, j_ind, edgeList)
        end
        t_map[i] = this_map
    end
    return( t_map )
end

# Initialize edgelist probability cache
function initialize_cache!(lcn::LatentChannelNetwork,
                           makeInverseMap::Bool = true)
    if makeInverseMap
        lcn.cache_map = makeTransposeMap(lcn.edgeList)
    end

    if (length(lcn.cache_map) == 0)
        error("map not initialized, set makeInverseMap == true")
    end

    nNodes = lcn.nNodes
    cache = Vector{Vector{Float64}}(undef, nNodes)
    for i in 1:nNodes
        these_edges = lcn.edgeList[i]
        these_probs = zeros(length(these_edges))
        for j in 1:length(these_edges)
            j_ind = these_edges[j]
            these_probs[j] = probEdge(i,j_ind, lcn)
        end
        cache[i] = these_probs
    end
    lcn.cache_probs = cache
end

# Find mapping
function findTransposeInds(i::Int64, j::Int64,
                           edgeList::Vector{Vector{Int64}})::Int64
    these_edges = edgeList[j]
    for ii in 1:length(these_edges)
        if these_edges[ii] == i
            return(ii)
        end
    end
    error("Lookup unsuccessful!")
end

# Get mapping
function getTransposeInds(i::Int64, j_ind::Int64,
                          lcn::LatentChannelNetwork)::Int64
    return(lcn.cache_map[i][j_ind])
end

# Insert updated probabilities
function updateCacheProb!(probs::Vector{Float64}, i::Int64,
                          lcn::LatentChannelNetwork)
    lcn.cache_probs[i] = probs
    i_revs = lcn.edgeList[i]
    for j_count in 1:length(probs)
        j_rev = getTransposeInds(i, j_count, lcn)
        lcn.cache_probs[i_revs[j_count]][j_rev] = probs[j_count]
    end
end

# Check that cache has been built
function checkCacheInit(lcn::LatentChannelNetwork)::Bool
    ans = (length(lcn.cache_probs) > 0)
    return(ans)
end



# Makes an latent channel network
function makeLatChan(edgeList::Array{Int64, 2},
                     dims::Int64 = 2)::LatentChannelNetwork

    preppedEdgeList = prepEdges(edgeList)
    # Getting number of nodes from edgeList
    nNodes = length(preppedEdgeList)
    # Random starting coordinates
    probs = rand(nNodes, dims)
    prob_means = vec( mean(probs, dims = 1) )
    empty_cache_probs = Vector{Vector{Float64}}(undef, 0)
    empty_map = Vector{Matrix{Int64}}(undef, 0)
    # Making Latent Hub Modeler
    ans = LatentChannelNetwork(nNodes, dims,
                               preppedEdgeList,
                               probs, prob_means,
                               empty_cache_probs, empty_map)

    initialize_cache!(ans)

    return(ans)
end

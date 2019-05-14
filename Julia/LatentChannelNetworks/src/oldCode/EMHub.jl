# include("edgeListUtils.jl")

using Statistics


# Latent Hub Class
mutable struct EMHub
    # Number of Nodes in Graph
    nNodes::Int64
    # Dimension of embedding
    dim::Int64
    # List of pairs
    edgeList::Array{Array{Int64, 1}}
    # Coordinates of Embedded Nodes
    hub_prob::Matrix{Float64}
    # Hub probability averages
    hub_prob_mean::Vector{Float64}

    # Cached table of edge probabilities for observed edges
    cache_probs::Vector{Vector{Float64}}
    # Reverse lookup tool so we can easily update i,j and j,i
    cache_map::Vector{Vector{Int64}}
end

# Makes an Graph Embedder with random parameter values
function makeLatHub(edgeList::Array{Int64, 2}, dims::Int64 = 2)::EMHub
    preppedEdgeList = prepEdges(edgeList)
    # Getting number of nodes from edgeList
    nNodes = length(preppedEdgeList)
    # Random starting coordinates
    probs = rand(nNodes, dims)
    prob_means = vec( mean(probs, dims = 1) )
    empty_cache_probs = Vector{Vector{Float64}}(undef, 0)
    empty_map = Vector{Matrix{Int64}}(undef, 0)
    # Making Latent Hub Modeler
    ans = EMHub(nNodes, dims, preppedEdgeList, probs, prob_means,
                empty_cache_probs, empty_map)
    return(ans)
end

# Compute probability of edge between nodes i and j
function probEdge(i::Int64, j::Int64, eh::EMHub)::Float64
    prob_no_edge = 1.0
    nDims = eh.dim
    probs = eh.hub_prob
    for ii in 1:nDims
        prob_no_edge = prob_no_edge * (1.0 - probs[i, ii] * probs[j,ii])
    end
    ans = 1.0 - prob_no_edge
    return(ans)
end

# For one step of EM algorithm, need to compute
# probability two nodes have edge, but exclude case in
# a certain hub connection is responsible for edge
function removeEdgeProbContribution(pik::Float64, pjk::Float64,
                                    edgeProb::Float64)::Float64
    probNoEdge = 1.0 - edgeProb
    probNoEdge_removePairCont = probNoEdge / (1.0 - pik * pjk)
    ans = 1.0 - probNoEdge_removePairCont
    return(ans)
end

function addEdgeProbContribution(pik::Float64, pjk::Float64,
    edgeProb::Float64)::Float64
    probNoEdge = 1.0 - edgeProb
    newProbNoEdge = probNoEdge * (1.0 - pik * pjk)
    ans = 1. - newProbNoEdge
    return(ans)
end

function updateEdgeProb(pik_new::Float64, pik_old::Float64,
    pjk::Float64, edgeProb::Float64)::Float64
    probNoEdge = 1.0 - edgeProb
    ans = 1.0 - probNoEdge / (1. - pik_old * pjk) * (1. - pik_new * pjk)
    return(ans)
end

# Extracts a boolean vector of whether node is connected to others
function getEdgeBools(i::Int64, eh::EMHub)::Vector{Bool}
    # Getting number of nodes
    nNodes = eh.nNodes
    # Creating vector of booleans
    ans = fill(false, nNodes)
    # Extracting edges for this node
    these_edges = eh.edgeList[i]
    # Setting nodes equal to true
    for j in 1:length(these_edges)
        ans[these_edges[j]] = true
    end
    return(ans)
end

# Compute Log-Likelihood contribution from single node
function computeNodeLLK(i::Int64, eh::EMHub)::Float64
    nNodes = eh.nNodes
    ans = 0.0
    these_edges = eh.edgeList[i]
    hasEdgeVec = getEdgeBools(i, eh)
    for j in 1:nNodes
        if i == j
            continue
        end
        this_edgeProb = probEdge(i,j,eh)
        if hasEdgeVec[j]
            ans = ans + log(this_edgeProb)
        else
            ans = ans + log(1.0 - this_edgeProb)
        end
    end
    return(ans)
end


# Computes Log Likelihood
function computeLLK(eh::EMHub)::Float64
    nNodes = eh.nNodes
    ans = 0.0
    for i in 1:nNodes
        ans = ans + computeNodeLLK(i, eh)
    end
    ans = ans / 2.0
    return(ans)
end


# Compute expected probability
# that node i has an edge into hub k going to
# node j, conditional on observing whether there is an edge between i and j
function expHubConnect(i::Int64, j::Int64, k::Int64,
                       hasEdge::Bool, eh::EMHub)::Float64
    p_ik = eh.hub_prob[i,k]
    p_jk = eh.hub_prob[j,k]
    if hasEdge
        edgeProb = probEdge(i,j,eh)
        edgeProb_contRemoved =
            removeEdgeProbContribution(p_ik, p_jk, edgeProb)

        top = p_ik * p_jk + p_ik * (1.0 - p_jk) * edgeProb_contRemoved
        ans = top / edgeProb
        return(ans)
    else
        ans = p_ik * (1.0 - p_jk)
        return(ans)
    end
end

# Updates a single edge probability
function em_oneProb!(i::Int64, k::Int64, eh::EMHub)
    hasEdgeVec = getEdgeBools(i, eh)
    prob_sum = 0.0
    nNodes = length(hasEdgeVec)
    for j in 1:nNodes
        if i == j
            continue
        end
        prob_sum = prob_sum + expHubConnect(i,j,k, hasEdgeVec[j], eh)
    end
    prob_mean = prob_sum / (nNodes - 1.0)
    eh.hub_prob[i,k] = prob_mean
end

# Updates single edge probability, but faster!
function em_oneProb_fast!(i::Int64, k::Int64, eh::EMHub)
    # Extracting information about this probability
    this_i_prob = eh.hub_prob[i,k]
    these_edges = eh.edgeList[i]
    nNodes = eh.nNodes
    # Computing "naive" sum (i.e., if no edges)
    prob_sum = nNodes * this_i_prob * ( 1. - eh.hub_prob_mean[k] )
    # Removing contribution from own node to sum
    prob_sum = prob_sum - this_i_prob * (1. - this_i_prob)
    # Readjusting sum for edges
    for j in these_edges
        prob_sum = prob_sum - this_i_prob * (1. - eh.hub_prob[j,k])
        # Substantial improvement could be made
        # by caching edge probabilities in step below
        prob_sum = prob_sum + expHubConnect(i,j,k, true, eh)
    end

    # Computing new probability
    prob_mean = prob_sum / (nNodes - 1)
    # Computing change in probability to
    # efficiently update new mean probabilities
    prob_diff = prob_mean - this_i_prob
    # Updating new probability
    eh.hub_prob[i,k] = prob_mean
    # Updating average probability
    eh.hub_prob_mean[k] = eh.hub_prob_mean[k] + prob_diff / nNodes
end


# EM update on all probabilities
function em_all!(eh::EMHub, fast = true, update_means = false)
    nDim = eh.dim
    nNodes = eh.nNodes
    if update_means
        eh.hub_prob_mean = vec(mean(eh.hub_prob, dims = 1) )
    end
    for i in 1:nNodes
        for k in 1:nDim
            if fast
                em_oneProb_fast!(i,k,eh)
            else
                em_oneProb!(i,k, eh)
            end
        end
    end
end

function em_all!(eh::EMHub, iters::Int, fast = true)
    for i in 1:iters
        em_all!(eh, fast)
    end
end



# Takes a MH step, with a uniform proposal step
function MH_unif_node!(i::Int64, eh::EMHub)
    old_llk = computeNodeLLK(i, eh)
    old_vals = eh.hub_prob[i,:]
    nVals = length(old_vals)
    new_vals = rand(nVals)
    eh.hub_prob[i,:] = new_vals
    new_llk = computeNodeLLK(i,eh)

    ratio = exp(new_llk - old_llk)
    if ratio > rand()
        return
    else
        eh.hub_prob[i,:] = old_vals
    end
end

function MH_unif_all!(eh::EMHub)
    for i in 1:eh.nNodes
        MH_unif_node!(i, eh)
    end
end

# Computes probability of hub connections, conditional on graph having an edge
function getConditionalHubProb(i::Int64, j::Int64, eh::EMHub)
    edgeProb = probEdge(i,j,eh)
    K = eh.dim
    ans = zeros(K)
    for k in 1:K
        prob_connect = eh.hub_prob[i,k] * eh.hub_prob[j,k]
        ans[k] = prob_connect / edgeProb
    end
    return(ans)
end

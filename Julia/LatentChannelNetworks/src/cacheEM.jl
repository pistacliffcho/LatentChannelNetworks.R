
# Updates single probability using caching EM algorithm
function em_oneProb_cache!(i::Int64, k::Int64,
                          edgeProbs::Vector{Float64},
                          lcn::LatentChannelNetwork,
                          ptol::Float64 = 10^-10)
    # Extracting information about this probability
    start_prob = lcn.pmat[i,k]
    if (start_prob < ptol)
        return
    end
    these_edges = lcn.edgeList[i]
    nNodes = lcn.nNodes
    # Computing "naive" sum (i.e., if no edges)
    prob_sum = nNodes * start_prob * ( 1. - lcn.pbar[k] )
    # Removing contribution from own node to sum
    prob_sum = prob_sum - start_prob * (1. - start_prob)
    # Readjusting sum for edges
    for j in 1:length(these_edges)
        j_ind = these_edges[j]
        j_prob = lcn.pmat[j_ind, k]
        prob_sum = prob_sum - start_prob * (1. - j_prob)
        # Substantial improvement could be made
        # by caching edge probabilities in step below
        prob_sum = prob_sum + expHubConnect_hasEdge(start_prob, j_prob, edgeProbs[j])
    end

    # Computing new probability
    new_prob = prob_sum / (nNodes - 1)
    # Computing change in probability to
    # efficiently update new mean probabilities
    prob_diff = new_prob - start_prob
    # Updating new probability
    lcn.pmat[i,k] = new_prob
    # Updating average probability
    lcn.pbar[k] = lcn.pbar[k] + prob_diff / nNodes

    for j in 1:length(these_edges)
        j_ind = these_edges[j]
        edgeProbs[j] = updateEdgeProb(new_prob, start_prob,
                                      lcn.pmat[j_ind, k],
                                      edgeProbs[j])
    end
end

function em_oneSet_cache!(i::Int64, lcn::LatentChannelNetwork)
    edgeProbs = get_cache_probs(i, lcn)
    K = lcn.dim
    for k in 1:K
        em_oneProb_cache!(i, k, edgeProbs, lcn)
    end
    set_cache_probs!(edgeProbs, i, lcn)
end

function em_cached_iter!(lcn::LatentChannelNetwork)
    if !checkCacheInit(lcn)
        error("Cache not initialized")
    end
    nNodes = length(lcn.edgeList)
    for i in 1:nNodes
        em_oneSet_cache!(i, lcn)
    end
end

function em_cached!(lcn::LatentChannelNetwork,
                    iters::Int64 = 10000,
                    tol::Float64 = 0.0001,
                    warn::Bool = true)::Dict
    m_p_diff = tol + 1
    for i in 1:iters
        p_old = copy(lcn.pmat)
        em_cached_iter!(lcn)
        p_diff = abs.( lcn.pmat .- p_old )
        m_p_diff = maximum(p_diff)
        if (m_p_diff < tol)
            ans = Dict("its" => i, "p_diff" => m_p_diff)
            return(ans)
        end
    end
    if warn
        println("Warning: maximum iterations reached")
    end
    ans = Dict("its" => iters, "p_diff" => m_p_diff)
    return(ans)
end


# Run EM from several starting points, chose best
function multi_em!(lcn::LatentChannelNetwork,
                   nStarts::Int64 = 10,
                   start_iters = 500,
                   finish_iters = 10000,
                   verbose = true)::Dict
    pmat_dim = size(lcn.pmat)
    max_llk = -Inf
    pmat_max = copy(lcn.pmat)
    for i in 1:nStarts
        rand_mat = rand(pmat_dim[1], pmat_dim[2])
        lcn.pmat = rand_mat
        initialize_cache!(lcn, false)
        res = em_cached!(lcn, start_iters, 0.00001, false)
        this_llk = computeLLK(lcn)
        if verbose
            println("This log likelihood = ", this_llk)
        end
        if (this_llk > max_llk)
            max_llk = this_llk
            pmat_max = copy(lcn.pmat)
        end
    end
    lcn.pmat = pmat_max
    initialize_cache!(lcn, false)
    res = em_cached!(lcn, finish_iters)
    if verbose
        println("Final LLK = ", computeLLK(lcn))
    end
    return(res)
end

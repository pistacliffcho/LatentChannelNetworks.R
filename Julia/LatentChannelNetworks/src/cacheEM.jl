# Updates single probability using EM algorithm
function update_prob(i::Int64, k::Int64,
                     lcn::LatentChannelNetwork,
                     ptol::Float64 = 10^-10)::Float64
    # Extracting information about this probability
    pik_old = lcn.pmat[i,k]
    if (pik_old < ptol)
        return(pik_old)
    end
    these_edges = lcn.edgeList[i]
    these_edge_probs = lcn.cache_probs[i]
    nNodes = lcn.nNodes
    # Computing "naive" sum (i.e., if no edges)
    prob_sum = nNodes * pik_old * ( 1. - lcn.pbar[k] )
    # Removing contribution from own node to sum
    prob_sum = prob_sum - pik_old * (1. - pik_old)
    # Readjusting sum for edges
    for j in 1:length(these_edges)
        j_ind = these_edges[j]

        pjk = lcn.pmat[j_ind, k]
        # Subtracting off contribution of not having an edge
        prob_sum = prob_sum - pik_old * (1. - pjk)
        # Grabbing cached edge prob
        pi_ij = these_edge_probs[j]
        # Adding contribution of having an edge
        prob_sum = prob_sum + edgeSumContribution(pik_old, pjk, pi_ij)
    end
    # Computing new probability
    pik_new = prob_sum / (nNodes - 1.)
    return(pik_new)
end



function em_oneProb_cache!(i::Int64, k::Int64,
                          lcn::LatentChannelNetwork,
                          ptol::Float64 = 10^-10)::Nothing
    pik_old = lcn.pmat[i,k]
    pik_new = update_prob(i, k, lcn, ptol)
    # Computing change in probability to
    # efficiently update new mean probabilities
    prob_diff = pik_new - pik_old
    # Updating new probability
    lcn.pmat[i,k] = pik_new
    # Updating average probability
    lcn.pbar[k] = lcn.pbar[k] + prob_diff / lcn.nNodes

    # Updating cached edge probabilities
    this_map = lcn.cache_map[i]
    these_edges = lcn.edgeList[i]
    these_edge_probs = lcn.cache_probs[i]
    for j in 1:length(these_edges)
        # Extracting pjk
        j_ind = these_edges[j]
        pjk = lcn.pmat[j_ind, k]
        # Extracting old edge prob
        old_edgeProb = these_edge_probs[j]
        # Updating edge prob with new values
        new_edgeProb = updateEdgeProb(pik_new, pik_old,
                                      pjk, old_edgeProb)
        # Updating cached prob
        these_edge_probs[j] = new_edgeProb
        # Getting transposed indices
        t_i = these_edges[j]
        t_j = this_map[j]
        # Updating transposed indices
        lcn.cache_probs[t_i][t_j] = new_edgeProb
    end
end

function em_cached_iter!(lcn::LatentChannelNetwork)
    if !checkCacheInit(lcn)
        error("Cache not initialized")
    end
    nNodes = lcn.nNodes
    K = lcn.dim
    for i in 1:nNodes
        for k in 1:K
            em_oneProb_cache!(i, k, lcn)
        end
    end
end

function getMaxDiff(m1::Matrix{Float64}, m2::Matrix{Float64})::Float64
    N = length(m1)
    max_diff = 0
    for i in 1:N
        diff = abs( m1[i] - m2[i] )
        max_diff = max(max_diff, diff)
    end
    return(max_diff)
end

function copyMatVals(matTo::Matrix{Float64}, matFrom::Matrix{Float64})
    N = length(matTo)
    for i in 1:N
        matTo[i] = matFrom[i]
    end
end

function em_cached!(lcn::LatentChannelNetwork,
                    iters::Int64 = 10000,
                    tol::Float64 = 0.0001,
                    warn::Bool = true)::Dict
    m_p_diff = tol + 1
    p_old = copy(lcn.pmat)
    for i in 1:iters
        copyMatVals(p_old, lcn.pmat)
        em_cached_iter!(lcn)
#        p_diff = abs.( lcn.pmat .- p_old )
        m_p_diff = getMaxDiff(lcn.pmat, p_old)
        if (m_p_diff < tol)
            ans = Dict("its" => i, "p_diff" => m_p_diff)
            return(ans)
        end
    end
    if warn
        println("Warning: maximum iterations reached")
    end
    ans = Dict("its" => iters, "err" => m_p_diff)
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
    initialize_cache!(lcn)
    res = em_cached!(lcn, finish_iters)
    if verbose
        println("Final LLK = ", computeLLK(lcn))
    end
    return(res)
end

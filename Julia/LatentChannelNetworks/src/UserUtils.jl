using Plots
using StatsBase

# Looking at difference in average hub strength by group
function hubDivergenceByLabel(this_label, label_vec,
                              lcn::LatentChannelNetwork)::Vector{Float64}
    is_label = this_label .== label_vec
    is_not_label = this_label .!= label_vec
    p_isLabel = copy(lcn.pmat[is_label,:])
    p_isNotLabel = copy(lcn.pmat[is_not_label,:])

    K = eh.dim
    ans = zeros(K)
    for k in 1:K
        phat_label = sum( p_isLabel[:,k] ) / length(p_isLabel)
        phat_notLabel = sum( p_isNotLabel[:,k] ) / length(p_isNotLabel)
        ans[k] = phat_label - phat_notLabel
    end
    return(ans)
end

# Makes a heatmap by label
function heatMap(lcn::LatentChannelNetwork,
                 label::Vector,
                 colorGradient = :pu_or,
                 divideColor = "red",
                 divideWidth = 0.5,
                 sortByHubprob = true,
                 drawDivide = true)
    new_order = sortperm(label)
    p_mat = copy(lcn.pmat[new_order,:])
    if sortByHubprob
        max_label = 1
        hub_prob_diffs = hubDivergenceByLabel(max_label, label, lcn)
        new_hub_order = sortperm(hub_prob_diffs, rev = true)
        p_mat = p_mat[:,new_hub_order]
    end
    p_mat = transpose(p_mat)
    ans = heatmap(p_mat, c = colorGradient)
    if drawDivide
        max_h = size(p_mat)[1]
        slab = sort(label)
        for i in 1:(length(slab) - 1)
            if (slab[i] != slab[i+1])
                plot!([i,i], [0.5, max_h + 0.5],
                      color = divideColor,
                      linewidth = divideWidth,
                      leg = false)
            end
        end
    end
    return(ans)
end


function relabelBySize(x::Vector)::Vector
    x_rank = denserank(x)
    cnt_map = countmap(x_rank)
    cnts = values(cnt_map)
    cnt_array = Matrix{Int64}(undef, length(cnts), 2)
    for (k,v) in cnt_map
        cnt_array[k,1] = k
        cnt_array[k,2] = v
    end
    vs = cnt_array[:,2]
    rnk = ordinalrank(vs, rev = true)
    ans = rnk[x_rank]
    return(ans)
end


# Compute the hub size
function hubSize(k::Int64, lcn::LatentChannelNetwork)::Float64
    ans = sum(lcn.pmat[:,k])
    return(ans)
end

# Compute all hub sizes
function hubSizes(lcn::LatentChannelNetwork)::Vector{Float64}
    ans = zeros(lcn.dim)
    for i in 1:lcn.dim
        ans[i] = hubSize(i, lcn)
    end
    return(ans)
end

# Compute conditional probability of connection through latent hub
function computeTheta(i::Int64,
                      j::Int64,
                      lcn::LatentChannelNetwork)::Vector{Float64}
    margEdgeProb = probEdge(i, j, lcn)
    K = eh.dim
    ans = zeros(K)
    for k in 1:K
        ans[k] = lcn.pmat[i,k] * lcn.pmat[j,k] / margEdgeProb
    end
    return(ans)
end


# Compute conditional expected number of connections through hub
function expHubConnections(i::Int64, lcn::LatentChannelNetwork)::Vector{Float64}
    ans = zeros(lcn.dim)
    K = length(ans)
    these_edges = lcn.edgeList[i]
    for e in these_edges
        ans = ans .+ computeTheta(i, e, lcn)
    end
    return(ans)
end

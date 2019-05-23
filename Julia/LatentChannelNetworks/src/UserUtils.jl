using Plots
using StatsBase

# Looking at difference in average hub strength by group
function hubDivergenceByLabel(this_label, label_vec,
                              lcn::LatentChannelNetwork)::Vector{Float64}
    is_label = this_label .== label_vec
    is_not_label = this_label .!= label_vec
    p_isLabel = copy(lcn.pmat[is_label,:])
    p_isNotLabel = copy(lcn.pmat[is_not_label,:])

    K = lcn.dim
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
                 label::Vector;
                 colorGradient = :pu_or,
                 divideColor = "red",
                 divideWidth = 0.5,
                 sortByHubprob = true,
                 drawDivide = true,
                 plot_size = (650, 450))
    new_order = sortperm(label)
    sorted_label = sort(label)

    # Deriving vectors for tick markers
    tick_name = []
    tick_location = []
    old_point = 0.0
    for i in 1:(length(sorted_label) - 1)
        if sorted_label[i] != sorted_label[i+1]
            push!(tick_name, sorted_label[i])
            push!(tick_location, (old_point + i) / 2.0)
            old_point = i + 1
        end
    end
    n = length(sorted_label)
    push!(tick_name, sorted_label[n])
    push!(tick_location, (old_point + n) / 2.0)

    p_mat = copy(lcn.pmat[new_order,:])
    if sortByHubprob
        max_label = 1
        hub_prob_diffs = hubDivergenceByLabel(max_label, label, lcn)
        new_hub_order = sortperm(hub_prob_diffs, rev = true)
        p_mat = p_mat[:,new_hub_order]
    end
    p_mat = transpose(p_mat)
    K = lcn.dim
    ans = heatmap(p_mat, c = colorGradient,
                  xticks = (tick_location, tick_name),
                  xrotation = 270,
                  yticks = collect(1:K),
                  size = plot_size)
    if drawDivide
        max_h = size(p_mat)[1]
        slab = sort(label)
        for i in 1:(length(slab) - 1)
            if (slab[i] != slab[i+1])
                x_loc = i - 0.5
                plot!([x_loc, x_loc], [0.5, max_h + 0.5],
                      color = divideColor,
                      linewidth = divideWidth,
                      leg = false)
            end
        end
    end
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
    K = lcn.dim
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

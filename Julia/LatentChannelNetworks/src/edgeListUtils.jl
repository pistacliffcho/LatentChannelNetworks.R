using DataFrames
using StatsBase


function dropLoops(edgeList::Array{Int, 2})
    isLoop = edgeList[:,1] .== edgeList[:,2]
    ans = edgeList[.!isLoop,:]
    return(ans)
end

# For efficiency, want to have edge pairs sorted by first
# node AND include copies of both directions, i.e.,
# if we have i,j pair, make a j,i pair and add to sorted list
function doubleAndSortEdges(edgeList::Array{Int, 2})
    col1 = edgeList[:,1]
    col2 = edgeList[:,2]
    nNodes = length(col1)
    doubledData = zeros(2 * nNodes, 2)
    # Adding data as entered
    doubledData[1:nNodes, :] = copy(edgeList)
    # Adding data with i,j switched
    doubledData[(nNodes+1):(2*nNodes),1] = copy(col2)
    doubledData[(nNodes+1):(2*nNodes),2] = copy(col1)

    col1_order = sortperm(doubledData[:,1])
    ans = doubledData[col1_order,:]
    return(ans)
end

# Takes in an Array of edges, prepares important info
function prepEdges(edgeList::Array{Int, 2})::Array{Array{Int64}}
    min_e = minimum(edgeList)
    if (min_e != 1)
        error("edgeList does not start with 1")
    end
    # Deleting loops
    edges_noLoops = dropLoops(edgeList)

    # Sorting edges
    sortedEdges = doubleAndSortEdges(edges_noLoops)
    # Dividing up by first node
    edge_df = DataFrame(n1 = sortedEdges[:,1], n2 = sortedEdges[:,2])
    split_list = groupby(edge_df, :n1)

    # Making edgeList, which is a vector of connected nodes
    ans = Array{Array}(undef, 0)
    maxNode = maximum(edgeList)

    emptyVec = Array{Int64}(undef, 0)

    currCount = 0
    nSplitData = length(split_list)

    for i in 1:maxNode
        currCount = currCount + 1
        if currCount > nSplitData
            push!(ans, emptyVec)
        else
            this_data  = split_list[currCount]
            this_node = this_data[1,1]
            if this_node != i
                currCount = currCount - 1
                push!(ans, emptyVec)
            else
                push!(ans, this_data[:,2])
            end
        end
    end
    # Removing duplicates
    for i in 1:length(ans)
        ans[i] = unique(ans[i])
    end
    return(ans)
end



function getConnectedNodes(node_id::Int, preppedEdgeList)::Array{Int, 1}
    ans = preppedEdgeList[node_id]
    return(ans)
end

# Randomly samples n nodes from a 1D array of nodes that are connected to node
function connectedNodeSample(n::Int,
                             connectedNodes::Array{Int, 1})::Array{Int, 1}
    ans = StatsBase.sample(connectedNodes, n, replace = false)
    return(ans)
end

# Randomly samples n nodes that are NOT connected to node
function unconnectedNodeSample(n::Int, max_node::Int,
                               connectedNodes::Array{Int, 1})::Array{Int, 1}
    nNodes = length( connectedNodes )
    # Sample cannot be a connected node OR this_node
    randFloats = rand(n) .* (max_node - nNodes - 1) .+ 1
    randInts = Int.(floor.(randFloats))
    for i in 1:n
        for j in 1:nNodes
            isGreater = randInts[i] >= connectedNodes[j]
            if isGreater
                randInts[i] = randInts[i] + 1
            end
        end
    end

    return(randInts)
end



# Simulates an edge list according to a Stochastic Block Matrix
# Randomly adds edge to any node that has none
function simSBM(nBlocks::Int = 5, nPerBlock::Int = 20, pIn::Real = 0.25, pOut::Real = 0.025)
    nvec1 = Int[]
    nvec2 = Int[]

    nTot = nBlocks * nPerBlock
    for i in 2:nTot
        for j in 1:(i-1)
            if mod(i, nBlocks) == mod(j, nBlocks)
                this_prob = pIn
            else
                this_prob = pOut
            end
            has_edge = rand(1)[1] < this_prob
            if has_edge
                push!(nvec1, i)
                push!(nvec2, j)
            end
        end
    end

    ans = hcat( nvec1, nvec2 )
    return(ans)
end

# Augment an SBM with a new group with new pIn/pOut
function augSBM(nObs::Int64, pIn::Float64,
                pOut::Float64,
                orgEdgeList::Array{Int, 2})
    ans = copy(orgEdgeList)
    maxOrgNode = maximum(ans)
    firstNew = maxOrgNode+1
    lastNew = maxOrgNode + nObs
    for i in firstNew:lastNew
        for j in 1:maxOrgNode
            if rand() < pOut
                new_pair = [i j;]
                ans = vcat(new_pair, ans)
            end
        end
        for j in firstNew:lastNew
            if i == j
                continue
            end
            if rand() < pIn
                new_pair = [i j;]
                ans = vcat(new_pair, ans)
            end
        end
    end
    return(ans)
end

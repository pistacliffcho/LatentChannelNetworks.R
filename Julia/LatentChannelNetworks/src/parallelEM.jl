# Serial implementation of simple EM algorithm
# This algorithm can be easily parallelized,
# but is expected to take slightly more
# iterations than em_cached
# Returns the error
function em_serial!(lcn::LatentChannelNetwork,
                  iters::Int64 = 10000,
                  tol::Float64 = 0.0001,
                  warn::Bool = true)::Dict
  iter = 0
  err = tol + 1
  for i in 1:iters
    err = serialEM_iter!(lcn)
    if (err < tol)
      ans = Dict("its" => i, "err" => err)
    end
  end
  if warn
    println("Warning: maximum iterations reached")
  end
  ans = Dict("its" =>iters, "err" => err)
  return(ans)
end

# Single iteration of serial EM algorithm
function serialEM_iter!(lcn::LatentChannelNetwork)::Float64
  pmat_new = copy(lcn.pmat)
  for i in 1:lcn.nNodes
    for k in 1:lcn.dim
      pmat_new[i,k] = update_prob(i,k,lcn)
    end
  end
  err = getMaxDiff(lcn.pmat, pmat_new)  
  lcn.pmat = pmat_new
  initialize_cache!(lcn)
  return(err)
end

# Initiailze edge probability cache in parallel
function parallel_init_cache!(lcn)
  Threads.@threads for i in 1:lcn.nNodes
    initialize_one_node!(i, lcn)
  end
  initialize_pbar!(lcn)
end


# Parallel implementation of simple EM algorithm
function em_parallel!(lcn::LatentChannelNetwork,
                  iters::Int64 = 10000,
                  tol::Float64 = 0.0001,
                  warn::Bool = true)::Dict
  iter = 0
  err = tol + 1
  for i in 1:iters
    err = parallel_EM_iter!(lcn)
    if (err < tol)
      ans = Dict("its" => i, "err" => err)
    end
  end
  if warn
    println("Warning: maximum iterations reached")
  end
  ans = Dict("its" =>iters, "err" => err)
  return(ans)
end

# Single iteration of parallel EM algorithm
function parallel_EM_iter!(lcn::LatentChannelNetwork)::Float64
  pmat_new = copy(lcn.pmat)
  Threads.@threads for i in 1:lcn.nNodes
    for k in 1:lcn.dim
      pmat_new[i,k] = update_prob(i,k,lcn)
    end
  end
  err = getMaxDiff(lcn.pmat, pmat_new)
  lcn.pmat = pmat_new
  parallel_init_cache!(lcn)
  return(err)
end

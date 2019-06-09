using Revise
using LatentChannelNetworks
using Test

println("Number of threads = ", Threads.nthreads())

nGrps = 10
nPerGrp = 100
pin = 0.25
pout = 0.05

sbm_edgeList = simSBM(nGrps, nPerGrp, pin, pout)
println("Number of edges = ", size(sbm_edgeList)[1],
        "\nNumber of Nodes = ", maximum(sbm_edgeList) ) 
lcn = makeLatChan(sbm_edgeList, nGrps)
starting_llk = computeLLK(lcn)
start = time()
@time res = em_cached!(lcn)
finish = time()
finishing_llk = computeLLK(lcn)
println("Finishing LLK for em_cached = ", finishing_llk)
println("Number of iterations = ", res["its"],
        " Number of seconds = ", round(finish - start))
println("Iterations per second = ", round(res["its"] / (finish - start)), "\n")

@test finishing_llk > starting_llk

# Rerunning with em_serial!
lcn = makeLatChan(sbm_edgeList, nGrps)
starting_llk = computeLLK(lcn)
start = time()
@time res = em_serial!(lcn)
finish = time()
finishing_llk = computeLLK(lcn)
println("Finishing LLK for em_serial = ", finishing_llk)
println("Number of iterations = ", res["its"],
        " Number of seconds = ", round(finish - start))
println("Iterations per second = ", round(res["its"] / (finish - start)), "\n")

@test finishing_llk > starting_llk

# Rerunning with em_parallel!
lcn = makeLatChan(sbm_edgeList, nGrps)
starting_llk = computeLLK(lcn)
start = time()
@time res = em_parallel!(lcn)
finish = time()
finishing_llk = computeLLK(lcn)
println("Finishing LLK for em_parallel = ", finishing_llk)
println("Number of iterations = ", res["its"],
        " Number of seconds = ", round(finish - start))
println("Iterations per second = ", round(res["its"] / (finish - start)), "\n")

@test finishing_llk > starting_llk

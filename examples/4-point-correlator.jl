using Revise
using ITensors
using ITensorCorrelators
using Random

using PyCall; ENV["KMP_DUPLICATE_LIB_OK"] = true
import PyPlot
const plt = PyPlot; 
plt.matplotlib.use("TkAgg"); ENV["MPLBACKEND"] = "TkAgg"; plt.pygui(true); plt.ion()

Random.seed!(1234)

function perm_ind(A, n, f)
    if n>1
        f = insertion_sort(A, n-1, f)
        val = A[n-1]
        j = n-2
        while (j>=1 && A[j] > val)
            f += 1
            A[j+1] = A[j]
            j = j-1
        end
        A[j+1] = val
    end
    return f
end

sizes = [100]
opers = [10,20,30,40]
t_MPO = zeros(length(opers))
t_opt = zeros(length(opers))

for (i, oper) in enumerate(opers)
    n = 100
    s = siteinds("S=1/2", n)
    psi = randomMPS(s, j -> isodd(j) ? "↑" : "↓"; linkdims=300)
    cor_ops = ("Z", "Z", "Z", "Z", "Z")
    #op_sites = [(1, 2, 3, 4), (1, 2, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5)]
    op_sites = Tuple{Vararg{Int}}[]
    #op_sites = NTuple{4,Int}[]
    for s in 1:oper
        aa = (rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)))
        if unique(aa) == [aa...]
            push!(op_sites, sort(aa))
        end
    end
    #op_sites = [(1, 1, 1, 1), (1, 1, 2, 2), (1, 2, 3, 3), (3, 3, 2, 1), (1, 2, 3, 4)]
    #op_sites = [(2,5,6,6), (2,5,6,7)]
    #op_sites = @show sort(op_sites)

    function correlator_MPO(psi, cor_ops, op_sites)
        sites = siteinds(psi)  
        C = Dict{NTuple{5,Int}, ComplexF64}()
        for l in op_sites 
            os = OpSum()
            os += cor_ops[1], l[1], cor_ops[2], l[2], cor_ops[3], l[3], cor_ops[4], l[4], cor_ops[5], l[5]
            corr = MPO(os, sites)
            C[l...] = inner(psi', corr, psi) #SC correlation function   
        end
        return C 
    end

    t_MPO[i] = @elapsed correlator_MPO(psi, cor_ops, op_sites)
    t_opt[i] = @elapsed correlator(psi, cor_ops, op_sites) #for the moment it only works with all four op on different sites
end

plt.figure(1)
plt.plot(opers, t_MPO, ".-", label = "MPO")
plt.plot(opers, t_opt, ".-", label = "Opt")
plt.legend()

    #round.(values(res) .- values(res_1), digits = 8)

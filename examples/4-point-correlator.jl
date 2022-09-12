using Revise
using ITensors
#include("../src/ITensorCorrelators.jl")
using .ITensorCorrelators
using Random

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

#=
cor_ops = ("X", "Y", "Z", "S+")
op_sites = [(5,6,7,4), (3,4,6,8), (7,6,2,4), (6,3,5,7), (5,6,7,8)]
indices = similar(op_sites); op_sorted=similar(op_sites)
for (i,op_site) in enumerate(op_sites)
    op_sorted[i] = tuple(sort([op_site...])...)
    indices[i] = tuple([findall(x->x==j,op_site)[1] for j in op_sorted[i]]...)
end

op_sites_dist = similar(unique(indices))
for i in 1:length(unique(indices))
    op_sites = @show [op_sorted[j] for j in findall(x->x==unique(indices)[i],indices)]
    cor_ops_ord = tuple([cor_ops[j] for j in unique(indices)[i]])
end

println(op_sorted)
println(indices)
f = perm_ind([op_sites...], length(op_sites)+1, 0)
println(f, op_sorted)
=#



n = 10
s = siteinds("S=1/2", n)
psi = randomMPS(s, j -> isodd(j) ? "↑" : "↓"; linkdims=300)
cor_ops = ("Z", "Z", "Z", "Z", "Z")
#op_sites = [(1, 2, 3, 4), (1, 2, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5)]
op_sites = Tuple{Vararg{Int}}[]
#op_sites = NTuple{4,Int}[]
for s in 1:100
    aa = (rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)))
    if unique(aa) == [aa...]
        push!(op_sites, sort(aa))
    end
end
#op_sites = [(1, 1, 1, 1), (1, 1, 2, 2), (1, 2, 3, 3), (3, 3, 2, 1), (1, 2, 3, 4)]
#op_sites = [(2,5,6,6), (2,5,6,7)]
op_sites = sort(op_sites)

function correlator_MPO(psi, cor_ops, op_sites)
    sites = siteinds(psi)  
    C = Dict{NTuple{5,Int}, ComplexF64}()
    for l in op_sites #multithreading on 16 threads
        os = OpSum()
        os += cor_ops[1], l[1], cor_ops[2], l[2], cor_ops[3], l[3], cor_ops[4], l[4], cor_ops[5], l[5]
        corr = MPO(os, sites)
        C[l...] = inner(psi', corr, psi) #SC correlation function   
    end
    return C 
end

res = @time correlator_MPO(psi, cor_ops, op_sites)
res_1 = @time correlator(psi, cor_ops, op_sites) #for the moment it only works with all four op on different sites

round.(values(res) .-values(res_1), digits = 8)

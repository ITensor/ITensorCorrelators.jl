using Revise
using ITensors
using ITensorCorrelators
using Random

Random.seed!(1234)

n = 70
s = siteinds("S=1/2", n)
psi = randomMPS(s, j -> isodd(j) ? "↑" : "↓"; linkdims=100)
cor_ops = ("Z", "Z", "Z", "Z")
#op_sites = [(1, 2, 3, 4), (1, 2, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5)]
op_sites = NTuple{4,Int}[]
for s in 1:2000
    aa = (rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)), rand(1:(n-1)))
    if unique(aa) == [aa...]
        push!(op_sites, sort(aa))
    end
end
op_sites = sort(op_sites)

function correlator_MPO(psi, cor_ops, op_sites)
    sites = siteinds(psi)  
    C = zeros(ComplexF64, length(psi), length(psi), length(psi), length(psi))
    for l in op_sites #multithreading on 16 threads
        os = OpSum()            
        os += cor_ops[1], l[1], cor_ops[2], l[2], cor_ops[3], l[3], cor_ops[4], l[4]
        corr = MPO(os, sites)
        C[l...] = inner(psi', corr, psi) #SC correlation function   
    end
    return C 
end

res = @time correlator_MPO(psi, cor_ops, op_sites)
res_1 = @time correlator(psi, cor_ops, op_sites) #for the moment it only works with all four op on different sites

res[findall(!=(0), res)] .- res_1[findall(!=(0), res_1)]
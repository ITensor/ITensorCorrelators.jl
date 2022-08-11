using Revise
using ITensors
using ITensorCorrelators
using Random

Random.seed!(1234)

n = 20
s = siteinds("S=1/2", n)
psi = randomMPS(s, j -> isodd(j) ? "↑" : "↓"; linkdims=40)
cor_ops = ("Z", "Z", "Z", "Z")
op_sites = [(1, 2, 3, 4), (1, 2, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5)]

function correlator_MPO(psi, cor_ops, op_sites)
    sites = siteinds(psi)  
    C = zeros(ComplexF64, length(psi), length(psi), length(psi), length(psi))
    for l in op_sites #multithreading on 16 threads
        os = OpSum()            
        os += cor_ops[1], l[1], cor_ops[2], l[2], cor_ops[3], l[3], cor_ops[4], l[4]
        println(os)           
        corr = MPO(os, sites)
        C[l...] = inner(psi', corr, psi) #SC correlation function   
    end
    return C 
end

res = correlator_MPO(psi, cor_ops, op_sites)
res_1 = correlator(psi, cor_ops, op_sites)

@show @time res[findall(!=(0), res)]
@show @time res_1[findall(!=(0), res_1)]
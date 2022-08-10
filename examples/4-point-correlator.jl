using ITensors
using ITensorCorrelators
using Random

Random.seed!(1234)

n = 8
s = siteinds("S=1/2", n)
ψ = randomMPS(s, j -> isodd(j) ? "↑" : "↓"; linkdims=4)
cor_ops = ("Z", "Z", "Z", "Z")
op_sites = [(1, 2, 3, 4), (1, 2, 4, 5), (1, 3, 4, 5), (2, 3, 4, 5)]
res = correlator(ψ, cor_ops, op_sites)

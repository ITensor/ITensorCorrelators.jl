using Revise
using ITensors
using ITensorCorrelators
using Random

using PyCall; ENV["KMP_DUPLICATE_LIB_OK"] = true
import PyPlot
const plt = PyPlot; 
plt.matplotlib.use("TkAgg"); ENV["MPLBACKEND"] = "TkAgg"; plt.pygui(true); plt.ion()

Random.seed!(2345)

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
opers = [20,30]
t_MPO = zeros(length(opers))
t_opt = zeros(length(opers))

#ITensors.enable_auto_fermion()



n = 10
s = siteinds("Electron", n; conserve_qns = true)
psi = randomMPS(s, j -> isodd(j) ? "Up" : "Dn", linkdims = 100)
#psi = productMPS(s, j -> isodd(j) ? "1" : "0")

function fh(n)
    ampo = OpSum()
    for i in 1:(n-1)
        ampo += 1, "Cdagup", i, "Cup", i+1
        ampo += 1, "Cdagup", i+1, "Cup", i
        ampo += 1, "Cdagdn", i, "Cdn", i+1
        ampo += 1, "Cdagdn", i+1, "Cdn", i

        ampo += 3, "Nup", i, "Ndn", i
    end
    return ampo
end

sweeps = Sweeps(10)
setmaxdim!(sweeps, 30)

os = fh(n)
H = MPO(os, s)
E, psi = dmrg(H, psi, sweeps)

cor_ops = ("Cdagup", "Cup")
op_sites = [(1, 2)]
#op_sites = Tuple{Vararg{Int}}[]
#op_sites = NTuple{4,Int}[]
#=
for s in 1:20
    aa = (rand(1:(n-1)), rand(1:(n-1)))
    if unique(aa) == [aa...]
        push!(op_sites, sort(aa))
    end
end
=#
@show op_sites

function correlator_MPO(psi, cor_ops, op_sites)
    sites = siteinds(psi)
    C = Dict{NTuple{2,Int}, ComplexF64}()
    for l in op_sites 
        os = OpSum()
        os += cor_ops[1], l[1], cor_ops[2], l[2]#, cor_ops[3], l[3], cor_ops[4], l[4]
        corr = MPO(os, sites)
        C[l...] = inner(psi', corr, psi) #SC correlation function   
    end
    return C 
end

t_MPO = @time correlator_MPO(psi, cor_ops, op_sites)
t_opt = @time correlator(psi, cor_ops, op_sites) #for the moment it only works with all four op on different sites
display(t_MPO)
display(t_opt)
t = round.(values(t_MPO) .- values(t_opt), digits = 8)


#==
plt.figure(1)
plt.plot(opers, t_MPO, ".-", label = "MPO")
plt.plot(opers, t_opt, ".-", label = "Opt")
plt.legend()
==#

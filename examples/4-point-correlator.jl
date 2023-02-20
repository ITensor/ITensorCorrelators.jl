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

#ITensors.enable_auto_fermion()


function square_lattice_nn(Nx::Int, Ny::Int; kwargs...)::Lattice
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Ny > 2)
    N = Nx * Ny
    Nbond = 2N - Ny + (yperiodic ? 0 : -Nx)
    Nbond_nn = 2(Nx-1)*(Ny - (yperiodic ? 0 : 1))
    latt = Lattice(undef, Nbond + Nbond_nn)
    b = 0
    for n in 1:N
      x = div(n - 1, Ny) + 1
      y = mod(n - 1, Ny) + 1

      #nearest neighbour
      if x < Nx
        latt[b += 1] = LatticeBond(n, n + Ny, x, y, x + 1, y, "1")
      end
      if Ny > 1
        if y < Ny
          latt[b += 1] = LatticeBond(n, n + 1, x, y, x, y + 1, "1")
        end
        if yperiodic && y == 1
          latt[b += 1] = LatticeBond(n, n + Ny - 1, x, y, x, y + Ny - 1, "1")
        end
      end

      #next nearest neighbour
      if x < Nx && Ny > 1
        if y < Ny
            latt[b += 1] = LatticeBond(n, n + Ny + 1, x, y, x + 1, y + 1, "2")
        end
        if y > 1
            latt[b += 1] = LatticeBond(n, n + Ny - 1, x, y, x + 1, y - 1, "2")
        end
        if yperiodic && y == Ny
            latt[b += 1] = LatticeBond(n, n + 1, x, Ny, x + 1, 1, "2")
        end
        if yperiodic && y == 1
            latt[b += 1] = LatticeBond(n, n + 2Ny - 1, x, 1, x + 1, Ny, "2")
        end
      end
    end
    return latt
end

function fh2(;t = 1, tp = 0.2, U = 10, Nx = 4, Ny = 2, α = 0.0)
    ampo = OpSum()  
    lattice = square_lattice_nn(Nx, Ny)
    for bond in lattice
        if bond.type == "1"
            pf = bond.y2 - bond.y1
            if abs(pf) > 1.1; pf = -sign(pf) end
            X = bond.x1 - (Nx)/2
            ampo += -t*exp(1im*2pi*α*X*pf), "Cdagup", bond.s1, "Cup", bond.s2 #nearest-neighbour hopping
            ampo += -t*exp(1im*2pi*α*X*pf), "Cdagdn", bond.s1, "Cdn", bond.s2
            ampo += -t*exp(-1im*2pi*α*X*pf), "Cdagup", bond.s2, "Cup", bond.s1
            ampo += -t*exp(-1im*2pi*α*X*pf), "Cdagdn", bond.s2, "Cdn", bond.s1
        end
        if bond.type == "2"
            pf = bond.y2-bond.y1
            if abs(pf)>1.1; pf = -sign(pf) end
            X = bond.x1 - (Nx)/2
            ampo += -tp*exp(1im*2pi*α*(X+0.5)*pf), "Cdagup", bond.s1, "Cup", bond.s2 #next-nearest-neighbour hopping
            ampo += -tp*exp(1im*2pi*α*(X+0.5)*pf), "Cdagdn", bond.s1, "Cdn", bond.s2
            ampo += -tp*exp(-1im*2pi*α*(X+0.5)*pf), "Cdagup", bond.s2, "Cup", bond.s1
            ampo += -tp*exp(-1im*2pi*α*(X+0.5)*pf), "Cdagdn", bond.s2, "Cdn", bond.s1
        end
    end
    for n in 1:(Nx*Ny)
        ampo += U, "Nup", n, "Ndn", n #density-density interaction
    end
    return ampo
end


function fh(n)
    ampo = OpSum()
    for i in 1:(n-1)
        ampo += -3, "Cdagup", i, "Cup", i+1
        ampo += -3, "Cdagup", i+1, "Cup", i
        #ampo += -1, "Cdagup", i, "Cdn", i+1
        #ampo += -1, "Cdagdn", i+1, "Cup", i
        ampo += -3, "Cdagdn", i, "Cdn", i+1
        ampo += -3, "Cdagdn", i+1, "Cdn", i

        ampo += 6, "Nup", i, "Ndn", i+1
    end
    return ampo
end

Nx = 5; Ny = 20
N = Nx * Ny
s = siteinds("Electron", N; conserve_qns = true)
psi = randomMPS(s, j -> isodd(j) ? "Up" : "Dn")
#psi = productMPS(s, j -> isodd(j) ? "1" : "0")


sweeps = Sweeps(10)
setmaxdim!(sweeps, 250)

#ss = fh2(; Nx = Nx, Ny = Ny)
ss = fh(N)
H = MPO(ss, s)
E, psi = dmrg(H, psi, sweeps)

cor_ops = ("Cdagup", "Cdagdn", "Cdn", "Cup")

op_sites = Tuple{Vararg{Int}}[]
for a in 1:(N-3)
    for b in (a+2):(N-1)
        push!(op_sites, (a, a+1, b, b+1))
    end
end

#op_sites = NTuple{4,Int}[]

#=
for s in 1:30
    a = rand(1:(n-1))
    b = rand(1:(n-1))
    aa = (a, a+1, b, b+1)
    if unique(aa) == [aa...]
        push!(op_sites, sort(aa))
    end
end 
=#

function correlator_MPO(psi, cor_ops, op_sites)
    sites = siteinds(psi)
    C = Dict{NTuple{4,Int}, ComplexF64}()
    for l in op_sites 
        os = OpSum()
        os += cor_ops[1], l[1], cor_ops[2], l[2], cor_ops[3], l[3], cor_ops[4], l[4]
        corr = MPO(os, sites)
        C[l...] = inner(psi', corr, psi) #SC correlation function   
    end
    return C 
end

t_MPO = @time correlator_MPO(psi, cor_ops, op_sites)
t_opt = @time correlator(psi, cor_ops, op_sites) #for the moment it only works with all four op on different sites
display(t_opt)
t = round.(values(t_MPO) .- values(t_opt), digits = 8)
display(t)
#==
plt.figure(1)
plt.plot(opers, t_MPO, ".-", label = "MPO")
plt.plot(opers, t_opt, ".-", label = "Opt")
plt.legend()
==#

#10  0.02    0.13   28
#20  0.15    2.16 153
#50 1.16  47 1128 
#100  4.61  418  4753

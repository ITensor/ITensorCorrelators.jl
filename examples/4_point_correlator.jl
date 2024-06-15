using ITensors, ITensorMPS
using ITensorCorrelators
using Random

Random.seed!(2345)

function ham(n; t, U)
  os = OpSum()
  for i in 1:(n - 1)
    os -= t, "Cdagup", i, "Cup", i + 1
    os -= t, "Cdagup", i + 1, "Cup", i
    os -= t, "Cdagdn", i, "Cdn", i + 1
    os -= t, "Cdagdn", i + 1, "Cdn", i
    os += U, "Nup", i, "Ndn", i + 1
  end
  return os
end

function correlator_mpo(psi, ops, op_sites)
  s = siteinds(psi)
  C = Dict{NTuple{4,Int},ComplexF64}()
  for l in op_sites
    os = OpSum()
    os += ops[1], l[1], ops[2], l[2], ops[3], l[3], ops[4], l[4]
    corr = MPO(os, s)
    C[l...] = inner(psi', corr, psi)
  end
  return C
end

function main(; N)
  s = siteinds("Electron", N; conserve_qns=true)
  psi = random_mps(s, j -> isodd(j) ? "Up" : "Dn")

  H = MPO(ham(N; t=1, U=2), s)
  E, psi = dmrg(H, psi; nsweeps=4, maxdim=50)

  ops = ("Cdagup", "Cdagdn", "Cdn", "Cup")

  sites = vec([(i, j, k, l) for i in 1:N, j in 1:N, k in 1:N, l in 1:N])

  correlators = @time correlator(psi, ops, sites)
  correlators_mpo = @time correlator_mpo(psi, ops, sites)
  return (; correlators, correlators_mpo)
end

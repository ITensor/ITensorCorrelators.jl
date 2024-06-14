using ITensors, ITensorMPS
using ITensorCorrelators

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

function main(N)
  sites = siteinds("Electron", N)
  ln = rand(1:50,length(sites)-1)
  psi = random_mps(sites, linkinds=ln) + im * random_mps(sites, linkinds=ln)

  # get all possible sets of indices
  indices = vec([tuple(i, j, k, l) for i in 1:N, j in 1:N, k in 1:N, l in 1:N])

  ops = ("Cdagup", "Cdagdn", "Cdn", "Cup")

  # using correlators
  correlators = @time correlator(psi, ops, indices)
  correlators_mpo = @time correlator_mpo(psi, ops, indices)

  return (; correlators, correlators_mpo)
end
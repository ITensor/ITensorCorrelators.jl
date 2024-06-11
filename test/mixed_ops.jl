
using ITensors, ITensorMPS
using Random
using Test
#using ITensorCorrelators
include("../src/ITensorCorrelators.jl")

function main()
  Random.seed!(0)

  sites = siteinds("Electron", 4)
  psi = randomMPS(sites, 3)

  op = ("Adagdn", "Adagup", "Cup", "Cdn")
  # cannot handle correlators with both fermionic and non-fermionic operators,
  # in particular, the sign is wrong when both fermionic and non-fermionic act on the same site
  inds = [(1, 3, 2, 1), (1, 2, 2, 3), (1, 2, 1, 3)]
  C = correlator(psi, op, inds)
  C2 = Dict{NTuple{4,Int},ComplexF64}()
  for idx in inds
    ampo = OpSum()
    ampo += op[1], idx[1], op[2], idx[2], op[3], idx[3], op[4], idx[4]
    mpo = MPO(ampo, sites)
    C2[idx] = inner(psi', mpo, psi)
  end
  return C == C2
end

function test()
  @test_broken main()
end

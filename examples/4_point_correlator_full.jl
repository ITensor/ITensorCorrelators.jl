using ITensors, ITensorMPS
using ITensorCorrelators

function main(N)
  sites = siteinds("Electron", N)
  psi = randomMPS(sites, 20) + im * randomMPS(sites, 20)

  # get all possible sets of indices
  indices = vec([tuple(i, j, k, l) for i in 1:N, j in 1:N, k in 1:N, l in 1:N])
  # using correlators
  C = correlator(psi, ("Cdagup", "Cdagdn", "Cdn", "Cup"), indices)
  # computing old-fashioned way
  for idx in indices
    ampo = OpSum()
    ampo += "Cdagup", idx[1], "Cdagdn", idx[2], "Cdn", idx[3], "Cup", idx[4]
    mpo = MPO(ampo, sites)
    # checking that correlators are the same
    @assert isapprox(inner(psi', mpo, psi), C[idx], rtol=1e-12, atol=1e-12)
  end
end

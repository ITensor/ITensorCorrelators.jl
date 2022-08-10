module ITensorCorrelators

using ITensors

# correlator(("A", "B", "C", "D"), [(1, 2, 3, 4), (1, 2, 4, 5), ...])
function correlator(
  psi::MPS,
  ops::NTuple{4,String},
  sites::Vector{NTuple{4,Int}},
)
  s = siteinds(psi)
  n = length(s)

  psidag = dag(sim(linkinds, psi))
  psi′ = psidag'

  # Sort the sites
  sites = sort(sites)

  # Check that each elements of `sites` is sorted
  @assert all(issorted, sites)

  #"Test"
  return sites
end

export correlator

end

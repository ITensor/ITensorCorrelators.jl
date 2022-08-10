module ITensorCorrelators

using ITensors

# correlator(("A", "B", "C", "D"), [(1, 2, 3, 4), (1, 2, 4, 5), ...])
function correlator(
  ops::NTuple{4,String},
  sites::Vector{NTuple{4,Int}},
  psi::MPS,
)
  s = siteinds(psi)
  n = length(s)
  sites = sort(sites)
  # Check that each elements of `sites` is sorted

end

export correlator

end

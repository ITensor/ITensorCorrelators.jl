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

  # Sort the sites
  sites = sort(sites)

  # Check that each elements of `sites` is sorted
  @assert all(issorted, sites)

  # for i, j, k, l
  # end

  is = unique(getindex.(sites, 1))

  @show is

  # for i in is
  i = is[1] # (i, j, k, l)
  sites_i = sites[findall(x -> x[1] == i, sites)]
  js = unique(getindex.(sites_i, 2))

  @show sites_i
  @show js

  o1_i = op(ops[1], s[i])
  L = apply(o1_i, psi[i]) * psidag[i]

  # for j in js
  j = js[1]

  for x in (i + 1):(j - 1)
    L = L * psi[x] * psidag[x]
  end

  return L
end

export correlator

end

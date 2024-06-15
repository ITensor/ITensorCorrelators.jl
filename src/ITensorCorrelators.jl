module ITensorCorrelators

using ITensors, ITensorMPS

include("correlator_bosonic_repeat.jl")
include("correlator_bosonic_recursive.jl")

function correlator(
  psi,
  cor_ops, #operator string to be computed, i.e. ("X","Y","Z")
  op_sites, #list of tuples with the sites where cor_ops are applied (does NOT have to be ordered)
)
  indices = similar(op_sites)
  op_sorted = similar(op_sites)

  for (i, op_site) in enumerate(op_sites)
    op_sorted[i] = tuple(sort([op_site...])...)
    indices[i] = tuple(sortperm([op_site...])...)
  end

  #C = Dict{Tuple{Vararg{Int64}}, ComplexF64}()
  #for i in 1:length(unique(indices)) #compute correlator separately for every possible order of sites
  i = 1
  ind_sites = unique(indices)[i]
  op_sites_ord = op_sites
  cor_ops_ord = cor_ops
  C = correlator_recursive_compact(psi, cor_ops_ord, op_sites_ord; indices=ind_sites)
  #end 
  return C
end

export correlator

end
#module ITensorCorrelators

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

  permutations = Dict{NTuple{length(op_sites[1]), Int64}, Vector{NTuple{length(op_sites[1]), Int64}}}()

  for (i, op_site) in enumerate(op_sites)
    op_sorted[i] = tuple(sort([op_site...])...)
    #indices[i] = tuple([findall(x -> x == j, op_site)[1] for j in op_sorted[i]]...) #store the order of the unordered sites strings
    indices[i] = tuple(sortperm([op_site...])...)
    op_sorted[i] in keys(permutations) ? push!(permutations[op_sorted[i]],indices[i]) : permutations[tuple(op_sorted[i]...)] = [indices[i]]
  end


  #C = Dict{Tuple{Vararg{Int64}}, ComplexF64}()
  #for i in 1:length(unique(indices)) #compute correlator separately for every possible order of sites
  i = 1
  ind_sites = unique(indices)[i]
  #ind_sites = indices[i]
  #op_sites_ord = [op_sorted[j] for j in findall(x -> x == ind_sites, indices)]
  op_sites_ord = op_sorted
  # if we only look at the first set of sites, we are not going to get all the orderings of the operators
  cor_ops_ord = tuple([cor_ops[j] for j in ind_sites]...)
  C = correlator_recursive_compact(psi, cor_ops_ord, op_sites_ord; indices=ind_sites,)
  #end 

  return C
end

export correlator

#end

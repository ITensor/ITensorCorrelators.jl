module ITensorCorrelators

using ITensors

include("correlator_bosonic.jl")
include("correlator_bosonic_repeat.jl")
include("correlator_fermionic.jl")
include("correlator_bosonic_recursive.jl")

function correlator(psi, cor_ops, op_sites)
    indices = similar(op_sites); op_sorted=similar(op_sites)
    for (i,op_site) in enumerate(op_sites)
        op_sorted[i] = tuple(sort([op_site...])...)
        indices[i] = tuple([findall(x->x==j,op_site)[1] for j in op_sorted[i]]...)
    end

    C = Dict{Tuple{Vararg{Int64}}, ComplexF64}()
    for i in 1:length(unique(indices))
        op_sites_ord = @show [op_sorted[j] for j in findall(x->x==unique(indices)[i],indices)]
        prinln(op_sites_ord)
        cor_ops_ord = tuple([cor_ops[j] for j in unique(indices)[i]])
        println(cor_ops_ord)
        push!(C,correlator_recursive_compact(psi, cor_ops_ord, op_sites_ord))
    end 
    return C
end

export correlator
export correlator_bosonic
export correlator_bosonic_repeat
export correlator_fermionic
export correlator_recursive
export correlator_recursive_compact


end

#TODO
#1) make it work for all hilbert spaces
#2) Add possibility of two operators on same site DONE
#3) Add variable numbers of loops
#4) add permutations 

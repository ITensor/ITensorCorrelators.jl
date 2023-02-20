module ITensorCorrelators

using ITensors

include("correlator_bosonic_repeat.jl")
include("correlator_bosonic_recursive.jl")

function fermi_to_bose(op)
    if op == "Cdag"
        return "Adag"
    elseif op == "C"
        return "A"
    elseif op == "Cdagup"
        return "Adagup"
    elseif op == "Cup"
        return "Aup"
    elseif op == "Cdagdn"
        return "Adagdn"
    elseif op == "Cdn"
        return "Adn"
    else 
        return op
    end
end

function correlator(psi,
    cor_ops, #operator string to be computed, i.e. ("X","Y","Z")
    op_sites #list of tuples witht the sites where cor_ops are applied (does NOT have to be ordered)
    )
    #cor_ops = map(a -> fermi_to_bose(a), cor_ops)
    indices = similar(op_sites); op_sorted = similar(op_sites)

    for (i,op_site) in enumerate(op_sites)
        op_sorted[i] = tuple(sort([op_site...])...)
        indices[i] = tuple([findall(x->x==j,op_site)[1] for j in op_sorted[i]]...) #store the order of the unordered sites strings
    end

    #C = Dict{Tuple{Vararg{Int64}}, ComplexF64}()
    #for i in 1:length(unique(indices)) #compute correlator separatedly for every possible order of sites
    i = 1
    ind_sites = unique(indices)[i]
    op_sites_ord = [op_sorted[j] for j in findall(x->x==ind_sites,indices)]
    cor_ops_ord = tuple([cor_ops[j] for j in ind_sites]...)
    C = correlator_recursive_compact(psi, cor_ops_ord, op_sites_ord; indices = ind_sites)
    #end 
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
#3) Add variable numbers of loops DONE
#4) add permutations DONE

using ITensors, ITensorMPS
using IterTools
using ProgressBars
using ITensorCorrelators

function main(N)
    operators = Dict([
        ("Qubit", ["X","Y","iY","Z","H","√NOT","Phase","π/8","Proj0","Proj1"]), # same as S=1/2
        ("S=1", ["Sz","Sz2","S+","S-","Sx","Sx2","iSy","Sy","Sy2"]),
        ("Boson", ["A","Adag","N"]),    # same as qudit
        ("Fermion", ["N","C","Cdag","F"]),
        ("Electron", ["Ntot","Nup","Ndn","Cup","Cdn","Cdagup","Cdagdn","Sz","Sx","S+","S-","F","Fup","Fdn"]),#,"Aup","Adagup","Adn","Adagdn"]),
        ("tJ", ["Ntot","Nup","Ndn","Cup","Cdn","Cdagup","Cdagdn","Sz","Sx","S+","S-","F","Fup","Fdn"])#,"Aup","Adagup","Adn","Adagdn"])
        ])
    for type in keys(operators) 
        println(type)
        sites = siteinds(type,N)
        psi = random_mps(sites,20) + im*random_mps(sites,20)
        
        # get all possible sets of indices
        indices = vec([tuple(i,j,k) for i=1:N,j=1:N,k=1:N])

        # construct all possible four point correlators
        ops = [operators[type] for i=1:3]
        for op in ProgressBar(IterTools.product(ops...))
            # check parity of fermionic operators
            ferms = [a in ["C","Cdag","Cup","Cdn","Cdagup","Cdagdn"] ? 1 : 0 for a in op]
            if sum(ferms) % 2 != 0
                continue
            else
                # using correlators
                C = correlator(psi,tuple(op...),indices)
                # computing old-fashioned way
                for idx in indices
                    ampo = OpSum()
                    ampo += op[1],idx[1],op[2],idx[2],op[3],idx[3]
                    mpo = MPO(ampo,sites)
                    # check the correlators are equivalent up to some numerical precision
                    @assert isapprox(inner(psi',mpo,psi),C[idx],rtol=1e-12,atol=1e-12) [op,idx]
                end
            end
        end
    end
end
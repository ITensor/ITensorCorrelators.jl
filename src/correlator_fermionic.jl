using ITensors

# correlator(("A", "B", "C", "D"), [(1, 2, 3, 4), (1, 2, 4, 5), ...])
function correlator_fermionic(
    psi::MPS,
    ops::NTuple{4,String},
    sites::Vector{NTuple{4,Int}},
    )
  
    sites = sort(sites) # Sort the sites
    @assert all(issorted, sites) # Check that each elements of `sites` is sorted
  
    C = Dict{NTuple{4,Int}, ComplexF64}()
    #C = zeros(ComplexF64, length(psi), length(psi), length(psi), length(psi)) #should we output a list or a NxNxNxN matrix?
  
    #psi = copy(psi0)    
    orthogonalize!(psi,1)
    psi_dag = prime(linkinds, dag(psi))
  
    is = unique(getindex.(sites, 1))
    for i in is #choose first site ("A")
        orthogonalize!(psi, i) #after orthogonalize weird things happen to the indices, have to do it before taking indices
        s = siteinds(psi) #this can be done before orth.
        ln = linkinds(psi)
        psi_dag = prime(linkinds, dag(psi)) #opposite MPS to contract
  
        sites_i = sites[findall(x -> x[1] == i, sites)] #checking the sites strings that start with i
        js = unique(getindex.(sites_i, 2)) #getting the corresponding js in position 2
  
        op_psi_i = apply(op(ops[1], s[i]), psi[i]) #apply first operator "A" in position i
  
        L_i = (i>1 ? delta(dag(ln[i-1])',ln[i-1]) : 1.) * op_psi_i * psi_dag[i]  #left system to conserve between contractions
        for str in (i+1):(js[1]-1) #contract between "A" and the first available "B"
          L_i = L_i * apply(op("F", s[str]), psi[str]) * psi_dag[str]
        end
  
        for (b, j) in enumerate(js) #choose second site ("B")
            sites_j = sites_i[findall(x -> x[2] == j, sites_i)] #checking the sites strings that has j in second position
            ks = unique(getindex.(sites_j, 3)) #getting the corresponding ks in position 3
  
            L_j = copy(L_i) #copy cached L_i
            op_psi_j =  apply(op(ops[2], s[j]),psi[j])  #apply second operator "B" in position j
            L_j = L_j * op_psi_j * psi_dag[j]  #generate left tensor to store
  
            for str in (j+1):(ks[1]-1) #contract between "B" and the first available site "C"
              L_j = L_j * psi[str] * psi_dag[str] 
            end
            for (c, k) in enumerate(ks) #choose third site ("C")
                sites_k = sites_j[findall(x -> x[3] == k, sites_j)] #checking the sites strings that has k in third position
                ls = unique(getindex.(sites_k, 4)) #getting the corresponding ls in position 4
  
                L_k = copy(L_j) #copy cached L_i
                op_psi_k = apply(op(ops[3], s[k]),psi[k])  #apply third operator "C" in position k
                L_k = L_k * op_psi_k * psi_dag[k]  #generate left tensor to store
  
                for str in (k+1):(ls[1]-1) #contract between "A" and the first available "B"
                  L_k = L_k * psi[str] * psi_dag[str] 
                end
                for (d, l) in enumerate(ls) #choose fourth site ("D")
  
                    L_l = copy(L_k) #copy cached L_k  
                    op_psi_l =  apply(op(ops[4],s[l]),psi[l]) #apply fourth operator "D" in position l
                    L_l = L_l * op_psi_l * psi_dag[l]  #generate left tensor to store
  
                    R = ((l)<length(psi) ? delta(dag(ln[l]),ln[l]') : 1.) #create right system
                    C[(i,j,k,l)] = inner(dag(L_l), R)
  
                    if (d<length(ls))
                      for str in (l):(ls[d+1]-1) #contract between "A" and the first available "B"
                        L_k = L_k * psi[str] * psi_dag[str] 
                      end 
                    end
                end
                if (c<length(ks))
                  for str in (k):(ks[c+1]-1) #contract between "A" and the first available "B"
                    L_j = L_j * psi[str] * psi_dag[str] 
                  end 
                end
            end
            if (b<length(js))
              for str in (j):(js[b+1]-1) #contract between "A" and the first available "B"
                L_i = L_i * psi[str] * psi_dag[str] 
              end 
            end
        end
    end
    return C
  end
  
@macroexpand @generated function mysum(A::Array{T,N}) where {T,N}
    quote
        s = zero(T)
        @nloops $N i A begin
            s += @nref $N A i
        end
        s
    end
end

module ITensorCorrelators

using ITensors

# correlator(("A", "B", "C", "D"), [(1, 2, 3, 4), (1, 2, 4, 5), ...])
function correlator(
  psi::MPS,
  ops::NTuple{4,String},
  sites::Vector{NTuple{4,Int}},
)
  sites = [(1, 2, 3, 4), (1, 2, 4, 5), (1, 3, 4, 5), (1, 3, 7, 8)]

  sites = sort(sites) # Sort the sites
  @assert all(issorted, sites) # Check that each elements of `sites` is sorted

  C = zeros(ComplexF64, length(lat), length(lat)) #should we output a list or a NxNxNxN matrix?

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
      L = (i>1 ? delta(dag(ln[i-1])',ln[i-1]) : 1.) * op_psi_i * psi_dag[i]  #left system to conserve between contractions

      for j in js #choose second site ("B")
          sites_j = sites[findall(x -> x[2] == j, sites)] #checking the sites strings that has j in second position
          ks = unique(getindex.(sites_j, 3)) #getting the corresponding ks in position 3

          L_j = copy(L) #copy cached L_i
          for str in (i+1):(j-1) #contract between "A" and "B"
            L_j = L_j * psi[str] * psi_dag[str] 
          end

          op_psi_j =  apply(op(ops[2], s[j]),psi[j])  #apply second operator "B" in position j

          L_j = L_j * op_psi_i1 * psi_dag[j]  #generate left tensor to store

          for k in ks #choose third site ("C")
              sites_k = sites[findall(x -> x[3] == k, sites)] #checking the sites strings that has k in third position
              ls = unique(getindex.(sites_k, 4)) #getting the corresponding ls in position 4
  
              L_k = copy(L_j) #copy cached L_i
              for str in (j+1):(k-1) #contract between "B" and "C"
                L_k = L_k * psi[str] * psi_dag[str] 
              end

              op_psi_k = apply(op(ops[3], s[k]),psi[k])  #apply third operator "C" in position k

              L_k = L_k * op_psi_k * psi_dag[k]  #generate left tensor to store
              for l in ls #choose fourth site ("D")

                  L_l = copy(L_k) #copy cached L_i
                  for str in (k+1):(l-1) #contract between "C" and "D"
                    L_l = L_l * psi[str] * psi_dag[str] 
                  end
  
                  op_psi_l =  apply(op(ops[4],s[l]),psi[l]) #apply fourth operator "D" in position l

                  L_l = L_l * op_psi_l * psi_dag[l]  #generate left tensor to store
                  R = ((l)<length(psi) ? delta(dag(ln[j+O]),ln[j+O]') : 1.) #create right system

                  C[i, j, k, l] = pf*inner(dag(L_l), R) #get matrix element
              end
              L_t = L_t * psi[j] * psi_dag[j]  #update left tensor (EFFICIENT PART)
          end
      end
  end


  return L
end

export correlator

end

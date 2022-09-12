using ITensors


function correlator_recursive_compact(
  psi::MPS,
  ops::Tuple{Vararg{String}},
  sites::Vector{Tuple{Vararg{Int}}},
  )

  sites = sort(sites) # Sort the sites
  N = length(sites[1])
  @assert all(issorted, sites) # Check that each elements of `sites` is sorted

  C = Dict{Tuple{Vararg{Int64}}, ComplexF64}()
  #C = zeros(ComplexF64, length(psi), length(psi), length(psi), length(psi)) #should we output a list or a NxNxNxN matrix?

  orthogonalize!(psi,1)
  psi_dag = prime(linkinds, dag(psi))
  op_inds = unique(getindex.(sites, 1))
  s = siteinds(psi) #this can be done before orth.
  ln = linkinds(psi)
  psi_dag = prime(linkinds, dag(psi)) #opposite MPS to contract
  L = ITensor(1.)
  counter = 1
  element = zeros(Int64, N)

  inner_module_compact(op_inds, sites, L, counter, element, N, ops, s, ln, psi, psi_dag, C)

  return C
end

function inner_module_compact(op_inds, sites_ind_prev, L_prev, counter, element, N, ops, s, ln, psi, psi_dag, C)
    for (a, op_ind) in enumerate(op_inds)
      element[counter] = op_ind
      if counter == 1
        orthogonalize!(psi, op_ind) #after orthogonalize weird things happen to the indices, have to do it before taking indices
        s = siteinds(psi) #this can be done before orth.
        ln = linkinds(psi)
        psi_dag = prime(linkinds, dag(psi)) #opposite MPS to contract
        L_prev = (op_ind>1 ? delta(dag(ln[op_ind-1])',ln[op_ind-1]) : 1.) 
      end

      L = copy(L_prev) #copy cached L_i
      op_psi =  apply(op(ops[counter], s[op_ind]),psi[op_ind])  #apply second operator "B" in position j
      L = L * op_psi * psi_dag[op_ind]  #generate left tensor to store
  
      if counter == N
        R = ((op_ind)<length(psi) ? delta(dag(ln[op_ind]),ln[op_ind]') : ITensor(1.)) #create right system
        C[tuple(element...)] = inner(dag(L), R)
      else
        sites_ind = sites_ind_prev[findall(x -> x[counter] == op_ind, sites_ind_prev)] #checking the sites strings that has j in second position
        op_inds_next = unique(getindex.(sites_ind, counter + 1)) #getting the corresponding ks in position 3  
  
        for str in (op_ind+1):(op_inds_next[1]-1) #contract between "B" and the first available site "C"
          L = L * psi[str] * psi_dag[str] 
        end
        inner_module_compact(op_inds_next, sites_ind, L, counter + 1, element, N, ops, s, ln, psi, psi_dag, C)
      end
  
      if (a<length(op_inds))
        for str in (op_ind):(op_inds[a+1]-1) #contract between "A" and the first available "B"
          L_prev = L_prev * psi[str] * psi_dag[str] 
        end 
      end
  
    end
  end
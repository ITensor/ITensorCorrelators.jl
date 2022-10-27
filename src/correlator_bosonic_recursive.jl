using ITensors

# correlator(("A", "B", "C", "D"), [(1, 2, 3, 4), (1, 2, 4, 5), ...])
function correlator_recursive(
    psi::MPS,
    ops::Tuple{Vararg{String}},
    sites::Vector{Tuple{Vararg{Int}}},
    )
  
    sites = sort(sites) # Sort the sites
    N = length(sites[1])
    @assert all(issorted, sites) # Check that each elements of `sites` is sorted
  
    C = Dict{Tuple{Vararg{Float64}}, ComplexF64}()
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
        counter = 1

        sites_i = sites[findall(x -> x[1] == i, sites)] #checking the sites strings that start with i
        js = unique(getindex.(sites_i, 2)) #getting the corresponding js in position 2
  
        op_psi_i = apply(op(ops[1], s[i]), psi[i]) #apply first operator "A" in position i
  
        L_i = (i>1 ? delta(dag(ln[i-1])',ln[i-1]) : 1.) * op_psi_i * psi_dag[i]  #left system to conserve between contractions
        for str in (i+1):(js[1]-1) #contract between "A" and the first available "B"
          L_i = L_i * psi[str] * psi_dag[str] 
        end

        add_operator(js, sites_i, L_i, counter + 1, N, ops, s, ln, psi, psi_dag)
        println("END")
      end
    return 0
  end

function add_operator(op_inds, sites_ind_prev, L_prev, counter, N, ops, s, ln, psi, psi_dag)
  for (a, op_ind) in enumerate(op_inds)
    L = copy(L_prev) #copy cached L_i
    op_psi =  apply(op(ops[counter], s[op_ind]),psi[op_ind])  #apply second operator "B" in position j
    L = L * op_psi * psi_dag[op_ind]  #generate left tensor to store

    if counter == N
      R = ((op_ind)<length(psi) ? delta(dag(ln[op_ind]),ln[op_ind]') : ITensor(1.)) #create right system
      @show inner(dag(L), R)
    else
      sites_ind = sites_ind_prev[findall(x -> x[counter] == op_ind, sites_ind_prev)] #checking the sites strings that has j in second position
      op_inds_next = unique(getindex.(sites_ind, counter + 1)) #getting the corresponding ks in position 3  

      for str in (op_ind+1):(op_inds_next[1]-1) #contract between "B" and the first available site "C"
        L = L * psi[str] * psi_dag[str] 
      end
      inner_module(op_inds_next, sites_ind, L, counter + 1, N, ops, s, ln, psi, psi_dag)
    end

    if (a<length(op_inds))
      for str in (op_ind):(op_inds[a+1]-1) #contract between "A" and the first available "B"
        L_prev = L_prev * psi[str] * psi_dag[str] 
      end 
    end

  end
end


function correlator_recursive_compact(
  psi, #::MPS,
  ops, #::Tuple{Vararg{String}},
  sites; #::Vector{Tuple{Vararg{Int}}},
  #indices = nothing
  )

  if indices === nothing
    indices = collect(1:length(ops))
  end

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

  add_operator(op_inds, sites, L, counter, element, N, ops, s, ln, psi, psi_dag, C, indices)

  return C
end

function add_operator(op_inds, sites_ind_prev, L_prev, counter, element, N, ops, s, ln, psi, psi_dag, C, indices)
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
        C[tuple([element[k] for k in indices]...)] = inner(dag(L), R)
      else
        sites_ind = sites_ind_prev[findall(x -> x[counter] == op_ind, sites_ind_prev)] #checking the sites strings that has j in second position
        op_inds_next = unique(getindex.(sites_ind, counter + 1)) #getting the corresponding ks in position 3  
  
        for str in (op_ind+1):(op_inds_next[1]-1) #contract between "B" and the first available site "C"
          L = L * psi[str] * psi_dag[str] 
        end
        add_operator(op_inds_next, sites_ind, L, counter + 1, element, N, ops, s, ln, psi, psi_dag, C, indices)
      end
  
      if (a<length(op_inds))
        for str in (op_ind):(op_inds[a+1]-1) #contract between "A" and the first available "B"
          L_prev = L_prev * psi[str] * psi_dag[str] 
        end 
      end
    end
  end
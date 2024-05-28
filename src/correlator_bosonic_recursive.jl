using TupleTools: TupleTools

# correlator(("A", "B", "C", "D"), [(1, 2, 3, 4), (1, 2, 4, 5), ...])
function correlator_recursive_compact(
  psi, #::MPS,
  ops, #::Tuple{Vararg{String}},
  sites; #::Vector{Tuple{Vararg{Int}}},
  indices=nothing,
)
  if indices === nothing #assumes all the sites are already properly ordered
    indices = collect(1:length(ops))
  end

  sites = sort(sites) # Sort the sites along y
  N = length(sites[1])
  @assert all(issorted, sites) # Check that each elements of `sites` is sorted

  C = Dict{Tuple{Vararg{Int64}},ComplexF64}() #initialize dictionary to store data
  #C = zeros(ComplexF64, length(psi), length(psi), length(psi), length(psi)) #should we output a list or a NxNxNxN matrix?

  orthogonalize!(psi, 1)
  psi_dag = prime(linkinds, dag(psi))
  @show sites
  op_inds = unique(getindex.(sites, 1))
  @show op_inds
  println("__________________________")
  s = siteinds(psi) #this can be done before orth.
  ln = linkinds(psi)
  psi_dag = prime(linkinds, dag(psi)) #opposite MPS to contract
  L = ITensor(1.0)
  counter = 1
  jw = 0 #keeps track of the number of fermionic operator to add a jordan-wigner term
  element = zeros(Int64, N)

  add_operator_fermi(
    op_inds, sites, L, counter, element, N, ops, s, ln, psi, psi_dag, C, indices, jw
  )
  return C
end

function add_operator(
  op_inds, sites_ind_prev, L_prev, counter, element, N, ops, s, ln, psi, psi_dag, C, indices
)
  for (a, op_ind) in enumerate(op_inds)
    element[counter] = op_ind
    if counter == 1
      orthogonalize!(psi, op_ind) #after orthogonalize weird things happen to the indices, have to do it before taking indices
      s = siteinds(psi) #this can be done before orth.
      ln = linkinds(psi)
      psi_dag = prime(linkinds, dag(psi)) #opposite MPS to contract
      L_prev = (op_ind > 1 ? delta(dag(ln[op_ind - 1])', ln[op_ind - 1]) : 1.0)
    end

    L = copy(L_prev) #copy cached L_i
    op_psi = apply(op(ops[counter], s[op_ind]), psi[op_ind])  #apply second operator "B" in position j
    L = L * op_psi * psi_dag[op_ind]  #generate left tensor to store

    if counter == N
      R = ((op_ind) < length(psi) ? delta(dag(ln[op_ind]), ln[op_ind]') : ITensor(1.0)) #create right system
      C[tuple([element[k] for k in [findall(x -> x == j, indices)[1] for j in TupleTools.sort(indices)]]...)] = inner(
        dag(L), R
      )
    else
      sites_ind = sites_ind_prev[findall(x -> x[counter] == op_ind, sites_ind_prev)] #checking the sites strings that has j in second position
      op_inds_next = unique(getindex.(sites_ind, counter + 1)) #getting the corresponding ks in position 3  

      for str in (op_ind + 1):(op_inds_next[1] - 1) #contract between "B" and the first available site "C"
        L = L * psi[str] * psi_dag[str]
      end
      add_operator(
        op_inds_next,
        sites_ind,
        L,
        counter + 1,
        element,
        N,
        ops,
        s,
        ln,
        psi,
        psi_dag,
        C,
        indices,
      )
    end
    if (a < length(op_inds))
      for str in (op_ind):(op_inds[a + 1] - 1) #contract between "A" and the first available "B"
        L_prev = L_prev * psi[str] * psi_dag[str]
      end
    end
  end
end

function add_operator_repeat(
  op_inds, sites_ind_prev, L_prev, counter, element, N, ops, s, ln, psi, psi_dag, C, indices
)
  for (a, op_ind) in enumerate(op_inds)
    element[counter] = op_ind
    if counter == 1
      orthogonalize!(psi, op_ind) #after orthogonalize weird things happen to the indices, have to do it before taking indices
      s = siteinds(psi) #this can be done before orth.
      ln = linkinds(psi)
      psi_dag = prime(linkinds, dag(psi)) #opposite MPS to contract
      L_prev = (op_ind > 1 ? delta(dag(ln[op_ind - 1])', ln[op_ind - 1]) : 1.0)
    end

    L = copy(L_prev) #copy cached L_i
    op_psi = apply(op(ops[counter], s[op_ind]), psi[op_ind])  #apply second operator "B" in position j
    L = L * op_psi * psi_dag[op_ind]  #generate left tensor to store

    if counter == N
      R = ((op_ind) < length(psi) ? delta(dag(ln[op_ind]), ln[op_ind]') : ITensor(1.0)) #create right system
      C[tuple([element[k] for k in [findall(x -> x == j, indices)[1] for j in TupleTools.sort(indices)]]...)] = inner(
        dag(L), R
      )
    else
      sites_ind = sites_ind_prev[findall(x -> x[counter] == op_ind, sites_ind_prev)] #checking the sites strings that has j in second position
      op_inds_next = unique(getindex.(sites_ind, counter + 1)) #getting the corresponding ks in position 3  

      for str in (op_ind + 1):(op_inds_next[1] - 1) #contract between "B" and the first available site "C"
        L = L * psi[str] * psi_dag[str]
      end
      add_operator(
        op_inds_next,
        sites_ind,
        L,
        counter + 1,
        element,
        N,
        ops,
        s,
        ln,
        psi,
        psi_dag,
        C,
        indices,
      )
    end
    if (a < length(op_inds))
      for str in (op_ind):(op_inds[a + 1] - 1) #contract between "A" and the first available "B"
        L_prev = L_prev * psi[str] * psi_dag[str]
      end
    end
  end
end

function add_operator_fermi(
  op_inds,
  sites_ind_prev,
  L_prev,
  counter,
  element,
  N,
  ops,
  s,
  ln,
  psi,
  psi_dag,
  C,
  indices,
  jw,
)
  for (a, op_ind) in enumerate(op_inds)
    element[counter] = op_ind 
    if counter == 1
      orthogonalize!(psi, op_ind) #after orthogonalize weird things happen to the indices, have to do it before taking indices
      s = siteinds(psi) #this can be done before orth.
      ln = linkinds(psi)
      psi_dag = prime(linkinds, dag(psi)) #opposite MPS to contract
      L_prev = (op_ind > 1 ? delta(dag(ln[op_ind - 1])', ln[op_ind - 1]) : 1.0) #initialize left environment
    end

    L = L_prev  #copy(L_prev) #copy cached environment after applying previous operator (this can be simplified)
    op_psi = psi[op_ind]

    if jw % 2 != 0
      op_psi = apply(op("F", s[op_ind]), op_psi) #apply jordan wigner string if needed
    end
    op_psi = apply(op(ops[counter], s[op_ind]), op_psi)  #apply operator in the spot #counter

    jw_next = 0
    if ops[counter] == "Cdagup" ||
      ops[counter] == "Cdag" ||
      ops[counter] == "Cup" ||
      ops[counter] == "C"
      jw_next = jw + 1
      op_psi = apply(op("F", s[op_ind]), op_psi) #track if a fermionic operator was applied
    elseif ops[counter] == "Cdagdn" || ops[counter] == "Cdn"
      jw_next = jw + 1
      op_psi = apply(op("F", s[op_ind]), op_psi) #for spin down operator we need a j-w term on-site
    end

    L = L * op_psi * psi_dag[op_ind]  #generate left environment to store

    if counter == N
      R = ((op_ind) < length(psi) ? delta(dag(ln[op_ind]), ln[op_ind]') : ITensor(1.0)) #create right system
      C[tuple([element[k] for k in [findall(x -> x == j, indices)[1] for j in TupleTools.sort(indices)]]...)] = inner(
        dag(L), R
      )
      #push!(C, tuple([element[k] for k in [findall(x->x==j,indices)[1] for j in sort(indices)]]...) => inner(dag(L), R))
      L = 0
    else
      @show sites_ind_prev
      sites_ind = sites_ind_prev[findall(x -> x[counter] == op_ind, sites_ind_prev)] #checking if there are more terms with the element #counter in operators to compute
      @show sites_ind
      op_inds_next = unique(getindex.(sites_ind, counter + 1)) #getting the sites counter+1 in the string 
      @show op_inds_next
      println("_____________________")
      for str in (op_ind + 1):(op_inds_next[1] - 1) #contract until the next operator to apply (with jw if required)
        if jw_next % 2 != 0
          L = L * apply(op("F", s[str]), psi[str]) * psi_dag[str]
        else
          L = L * psi[str] * psi_dag[str]
        end
      end
      add_operator_fermi(
        op_inds_next,
        sites_ind,
        L,
        counter + 1,
        element,
        N,
        ops,
        s,
        ln,
        psi,
        psi_dag,
        C,
        indices,
        jw_next,
      )
    end
    if (a < length(op_inds))
      for str in (op_ind):(op_inds[a + 1] - 1)
        if jw % 2 != 0
          L_prev = L_prev * apply(op("F", s[str]), psi[str]) * psi_dag[str]
        else
          L_prev = L_prev * psi[str] * psi_dag[str]
        end
      end
    end
  end
end

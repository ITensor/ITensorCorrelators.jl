using Combinatorics: parity

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

  C = Dict{Tuple{Vararg{Int64}},ComplexF64}() #initialize dictionary to store data
  #C = zeros(ComplexF64, length(psi), length(psi), length(psi), length(psi)) #should we output a list or a NxNxNxN matrix?

  orthogonalize!(psi, 1)
  psi_dag = prime(linkinds, dag(psi))

  # computes the first index, the number of repeats, and permutation of the operators
  inds_ord = [sort([sites[idx]...])[1] for idx in 1:length(sites)]
  repeats = [count(==(sort([sites[idx]...])[1]), sites[idx]) - 1 for idx in 1:length(sites)]
  perms = [sortperm([sites[idx]...])[1:(repeats[idx] + 1)] for idx in 1:length(sites)]

  op_inds = unique([(inds_ord[idx], repeats[idx], perms[idx]) for idx in 1:length(sites)])

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
    repeat = op_ind[2]  # counts how many times op_ind is repeated. repeat=0 means op_ind occurs once.
    perm_ind = op_ind[3]  # determines what operator acts on the site index
    op_ind = op_ind[1]  # the next site(s)

    element[counter:(counter + repeat)] .= op_ind
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

    for i in 0:repeat
      op_psi = apply(op(ops[perm_ind[counter + repeat - i]], s[op_ind]), op_psi)
    end

    jw_next = jw

    for i in 0:repeat
      if has_fermion_string(ops[perm_ind[counter + repeat - i]], s[op_ind])
        jw_next = jw_next + 1
        op_psi = apply(op("F", s[op_ind]), op_psi) #track if a fermionic operator was applied
      end
    end

    L = L * op_psi * psi_dag[op_ind]  #generate left environment to store

    if counter + repeat == N
      R = ((op_ind) < length(psi) ? delta(dag(ln[op_ind]), ln[op_ind]') : ITensor(1.0)) #create right system

      # re-arranging the sorted element list using permutations
      perm_elem = element[sortperm(perm_ind)]

      # checking for fermion operators and keeping track of anti-commutations
      ferm_sites =
        Int.(perm_elem[findall(x -> has_fermion_string(x, s[op_ind]), ops)])
      par = 1 - 2 * parity(sortperm(ferm_sites))

      C[tuple(perm_elem...)] = par * inner(dag(L), R)
      L = 0
    else

      # The filtering of sites in next iteration could probably be sped up

      # sets of sites consistent with the current site
      sites_ind = sites_ind_prev[findall(
        x -> sort([x...])[counter + repeat] == op_ind, sites_ind_prev
      )]

      # making sure the next site is not the same as the previous one, since the repeated indices are already taken care of
      deleteat!(
        sites_ind, findall(x -> sort([x...])[counter + repeat + 1] == op_ind, sites_ind)
      )

      # the number of repeats of the site in the next iteration
      repeat_next = [
        count(==(sort([sites_ind[idx]...])[counter + repeat + 1]), sites_ind[idx]) - 1 for
        idx in 1:length(sites_ind)
      ]

      # get the next sites and permutations 
      inds_ord = [
        sort([sites_ind[idx]...])[counter + repeat + 1] for idx in 1:length(sites_ind)
      ]
      perms = [
        sortperm([sites_ind[idx]...])[1:(counter + repeat + repeat_next[idx] + 1)] for
        idx in 1:length(sites_ind)
      ]

      # gather next site ind, the number of repeats, and permutation array into one vector
      op_inds_next = unique([
        (inds_ord[idx], repeat_next[idx], perms[idx]) for idx in 1:length(sites_ind)
      ])

      # making sure the next iteration has the same permutation vector up until this point
      op_inds_next = op_inds_next[findall(
        x -> x[3][1:length(perm_ind)] == perm_ind, op_inds_next
      )]

      for str in (op_ind + 1):(op_inds_next[1][1] - 1) #contract until the next operator to apply (with jw if required)
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
        counter + repeat + 1,
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
      for str in (op_ind):(op_inds[a + 1][1] - 1)
        if jw % 2 != 0
          L_prev = L_prev * apply(op("F", s[str]), psi[str]) * psi_dag[str]
        else
          L_prev = L_prev * psi[str] * psi_dag[str]
        end
      end
    end
  end
end
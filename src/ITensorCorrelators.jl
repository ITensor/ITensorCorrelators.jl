module ITensorCorrelators

using ITensors

include("correlator_bosonic.jl")
include("correlator_bosonic_repeat.jl")
include("correlator_fermionic.jl")
include("correlator_bosonic_recursive.jl")

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

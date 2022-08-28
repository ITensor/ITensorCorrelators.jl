module ITensorCorrelators

using ITensors

include("correlator_bosonic.jl")
include("correlator_bosonic_repeat.jl")
include("correlator_fermionic.jl")

export correlator_bosonic
export correlator_bosonic_repeat
export correlator_fermionic

end

#TODO
#1) make it work for all hilbert spaces
#2) Add possibility of two operators on same site DONE
#3) Add variable numbers of loops
#4) add permutations 

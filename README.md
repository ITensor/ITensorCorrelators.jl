# ITensorCorrelators

This is a package for efficiently measuring n-point correlators of matrix product states (MPS), built on top of the [ITensors.jl](https://github.com/ITensor/ITensors.jl) library.

## Installation

This package is currently not registered. You can install it with:
```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/ITensor/ITensorCorrelators.jl.git")
```

## Usage

As an example, to compute the correlator `<Sz_i Sz_j Sz_k Sz_l>` on some set of sites `i, j, k, l`, you can use the `correlator` function:
```julia
using ITensors
using ITensorCorrelators

s = siteinds("S=1/2", 10)
psi = randomMPS(s; linkdims=10)

c = correlator(psi, ("Sz", "Sz", "Sz", "Sz"), [(1, 3, 4, 5), (2, 3, 4, 5), (3, 4, 5, 10)])
c[(2, 3, 4, 5)]
```
This outputs a dictionary mapping the set of sites `i, j, k, l` to the resulting correlator `<Sz_i Sz_j Sz_k Sz_l>`, so in the last line above we are accessing the result for the correlator `<Sz_2 Sz_3 Sz_4 Sz_5>`.

Note that currently, the package does not handle the case of repeated sites properly (i.e. it doesn't handle inputs `(2, 3, 3, 5)` properly). However, it should handle fermionic operators like `Cdag` and `C` properly.

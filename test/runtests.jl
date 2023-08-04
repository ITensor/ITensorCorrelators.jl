using ITensorCorrelators
using Test

@testset "ITensorCorrelators.jl" begin
  include(joinpath(pkgdir(ITensorCorrelators), "examples", "4_point_correlator.jl"))
  (; correlators, correlators_mpo) = main(; N=20)
  op_sites = keys(correlators)
  vals = [correlators[op_site] for op_site in op_sites]
  vals_mpo = [correlators_mpo[op_site] for op_site in op_sites]
  @test vals ≈ vals_mpo rtol = 1e-12
end

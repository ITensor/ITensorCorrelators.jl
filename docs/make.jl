using ITensorCorrelators
using Documenter

DocMeta.setdocmeta!(
  ITensorCorrelators, :DocTestSetup, :(using ITensorCorrelators); recursive=true
)

makedocs(;
  modules=[ITensorCorrelators],
  authors="Matthew Fishman <mfishman@flatironinstitute.org> and contributors",
  repo="https://github.com/mtfishman/ITensorCorrelators.jl/blob/{commit}{path}#{line}",
  sitename="ITensorCorrelators.jl",
  format=Documenter.HTML(;
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://mtfishman.github.io/ITensorCorrelators.jl",
    edit_link="main",
    assets=String[],
  ),
  pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/mtfishman/ITensorCorrelators.jl", devbranch="main")

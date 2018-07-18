using Documenter, BEDFiles

makedocs(
    format = :html,
    sitename = "BEDFiles"
)

deploydocs(
    repo   = "github.com/dmbates/BEDFiles.jl.git",
    target = "build",
    julia  = "nightly",
    deps   = nothing,
    make   = nothing
)


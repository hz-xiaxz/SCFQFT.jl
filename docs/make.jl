using SCFQFT
using Documenter

DocMeta.setdocmeta!(SCFQFT, :DocTestSetup, :(using SCFQFT); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers

makedocs(;
    modules = [SCFQFT],
    authors = "hz-xiaxz <1806656034@qq.com> and contributors",
    repo = "https://github.com/hz-xiaxz/SCFQFT.jl/blob/{commit}{path}#{line}",
    sitename = "SCFQFT.jl",
    format = Documenter.HTML(; canonical = "https://hz-xiaxz.github.io/SCFQFT.jl"),
    pages = [
        "index.md"
        [
            file for file in readdir(joinpath(@__DIR__, "src")) if
            file != "index.md" && splitext(file)[2] == ".md"
        ]
    ],
)

deploydocs(; repo = "github.com/hz-xiaxz/SCFQFT.jl")

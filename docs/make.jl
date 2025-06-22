using PhasedArrays
using Documenter

DocMeta.setdocmeta!(PhasedArrays, :DocTestSetup, :(using PhasedArrays); recursive=true)

makedocs(;
    modules=[PhasedArrays],
    authors="Iaroslav Shilinkov <yashilinkov@gmail.com>",
    sitename="PhasedArrays.jl",
    format=Documenter.HTML(;
        canonical="https://Yashilinkov.github.io/PhasedArrays.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Yashilinkov/PhasedArrays.jl",
    devbranch="main",
)

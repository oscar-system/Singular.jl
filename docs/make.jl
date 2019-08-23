using Documenter, Singular

makedocs(
         format   = Documenter.HTML(),
         sitename = "Singular.jl",
         pages    = [
             "index.md",
             "Coefficient rings" => [
                  "integer.md",
                  "rational.md",
                  "modn.md",
                  "modp.md",
                  "GF.md",
                  "nemo.md"
             ],
            "polynomial.md",
            "ideal.md",
            "Modules" => [
                  "module.md",
                  "vector.md"
            ],
            "alghom.md",
            "resolution.md",
            "matrix.md"
         ]
)

deploydocs(
   repo   = "github.com/wbhart/Singular.jl.git",
   target = "build",
   deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
   make   = nothing
)


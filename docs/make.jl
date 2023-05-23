using Documenter
using Singular

makedocs(
         format = Documenter.HTML(),
         sitename = "Singular.jl",
         modules = [Singular],
         clean = true,
         doctest = false,
         pages    = [
             "index.md",
             "Coefficient rings" => [
                  "integer.md",
                  "rational.md",
                  "modn.md",
                  "modp.md",
                  "transExt.md",
                  "GF.md",
                  "nemo.md"
             ],
            "polynomial.md",
            "noncommutative.md",
            "ideal.md",
            "qring.md",
            "Modules" => [
                  "module.md",
                  "vector.md"
            ],
            "alghom.md",
            "resolution.md",
            "matrix.md",
            "caller.md"
         ]
)

deploydocs(
   repo   = "github.com/oscar-system/Singular.jl.git",
   target = "build",
   deps = nothing,
   make   = nothing,
)


using Documenter, Singular

makedocs(
         format   = :html,
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
            "polynomial.md"
         ]
)

deploydocs(
   julia = "nightly",
   repo   = "github.com/wbhart/Singular.jl.git",
   target = "build",
   deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
   make   = nothing
)


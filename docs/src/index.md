# Getting Started

Singular.jl is a Julia interface to the Singular computer algebra system. It was
written by Oleksandr Motsak, William Hart and other contributors, and is maintained by
William Hart, Hans Schoenemann and Andreas Steenpas. It is part of the Oscar project which is
funded by the DFG TRR 195 program.

- [https://www.singular.uni-kl.de/](https://www.singular.uni-kl.de/) (Singular website)
- [https://github.com/oscar-system/Singular.jl](https://github.com/oscar-system/Singular.jl) (Singular.jl source code)
- [http://oscar-system.github.io/Singular.jl/](http://oscar-system.github.io/Singular.jl/) (Singular.jl online documentation)

The features of Singular so far include:

  - Singular integers, rationals Z/nZ, Z/pZ, Galois fields
  - Multivariate polynomials
  - Ideals over polynomial rings
  - Free modules over polynomial rings and submodules given by a finite generating set
  - Groebner basis over a field
  - Free/minimal resolutions
  - Syzygy modules
  - Nemo.jl rings can be used as coefficient rings

## Installation

To use Singular.jl we require Julia 1.0 or higher. Please see
[http://julialang.org/downloads](http://julialang.org/downloads/) for instructions on
how to obtain julia for your system.

At the Julia prompt simply type

```
julia> using Pkg
julia> Pkg.clone("https://github.com/oscar-system/Singular.jl")
julia> Pkg.build("Singular")
```

Here is an example of using Singular.jl

```julia
julia> using Singular

julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(x,y),(dp(2),C), Singular.spoly{Singular.n_Q}[x, y])

julia> I = Ideal(R, x^2 + 1, x*y + 1)
Singular Ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x^2+1, x*y+1)

julia> G = std(I)
Singular Ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x-y, y^2+1)

julia> Z = syz(G)
Singular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:
y^2*gen(1)-x*gen(2)+y*gen(2)+gen(1)

julia> F = fres(G, 0)
Singular Resolution:
R^1 <- R^2 <- R^1

julia> F[1]
Singular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:
x-y
y^2+1
```


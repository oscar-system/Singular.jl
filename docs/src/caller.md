```@meta
CurrentModule = Singular
```

# Library Proceedures

Singular.jl provides limited access to the functionality of the Singular
interpreter and its associated standard library procedures. The following
features of the Singular language are not supported:

 - The singular language allows arbitrary attributes to be attached to each
   object. These are not supported because Singular.jl works with the raw
   kernel types ring, poly, ideal, etc. that do not have attributes.

 - Maps between polynomial rings cannot be effectively communicated between the
   Singular language and Singular.jl because the Singular language tracks the
   preimage ring by name only and the rings in Singular.jl do not have names.

 - Very few library proceedures work with [Nemo rings and fields](@ref)

## Calling a Library Procedure

In general, if we have a procedure `sort` in `general.lib`, then the
corresponding function in `Singular.jl` is called `Singular.LibGeneral.sort`.
The full list of libraries included can be viewed by typing `Singular.Lib` at
the REPL and double pressing the tab key.

One issue that comes up in calling library proceedures is the implicit argument
`basering` that all procedures receive in the Singular language. `Singular.jl`
tries to infer the base ring from the arguments provided to the function. When
this fails or is simply not possible, the user can always provide a base ring
by passing it in as the first argument to the `Singular.jl` function. Note that
if the first argument to the `Singular.jl` version of a library procedure, then
this is automatically assumed to be the base ring. Hence, if a procedure in the
Singular language takes a ring as a first argument, you will have to pass that
ring in as the second argument after specifying the base ring in the first
argument.

**Examples**

This example illustrates passing Singular lists and providing the base ring.

```julia
julia> r0, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"], ordering=ordering_lp());

julia> Singular.LibGeneral.sort([x, y])
ERROR: `intvec` may be passed in as Vector{Int}. All other vectors (`list` in Singular) must be passed in as Vector{Any} along with an explicit base ring in the first argument

julia> Singular.LibGeneral.sort(r0, Any[x, y])
2-element Vector{Vector{T} where T}:
 spoly{n_Q}[y, x]
 [2, 1]
```

This example illustrates the base ring inference:

```julia
julia> AA, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"]);

julia> D = zero_matrix(AA, 4, 4);

julia> D[1,2] = -z; D[1,3] = 2*x; D[2,3] = -2*y;

julia> A, (x, y, z, t) = GAlgebra(AA, 1, D);

julia> Singular.LibNctools.isCentral(x)   # base ring A is infered from x
0

julia> Singular.LibCentral.center(A, 3)   # base ring cannot be infered from the plain Int 3
Singular ideal over Singular G-Algebra (QQ),(x,y,z,t),(dp(4),C) with generators (t, 4*x*y + z^2 - 2*z)
```


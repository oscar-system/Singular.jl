```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Interpreter Functionality

Singular.jl provides limited access to the functionality of the Singular
interpreter and its associated standard library procedures. The following
features of the Singular language are not supported:

 - The singular language allows arbitrary attributes to be attached to each
   object. These are not supported because Singular.jl works with the raw
   kernel types ring, poly, ideal, etc. that do not have attributes.

 - Maps between polynomial rings cannot be effectively communicated between the
   Singular language and Singular.jl because the Singular language tracks the
   preimage ring by name only and the rings in Singular.jl do not have names.

 - Very few library procedures work with [Nemo rings and fields](@ref).

## Calling a Library Procedure

In general, if we have a procedure `sort` in `general.lib`, then the
corresponding function in `Singular.jl` is called `Singular.LibGeneral.sort`.
The full list of libraries included can be viewed by typing `Singular.Lib` at
the REPL and double pressing the tab key.

One issue that comes up in calling library procedures is the implicit argument
`basering` that all procedures receive in the Singular language. `Singular.jl`
tries to infer the base ring from the arguments provided to the function. When
this fails or is simply not possible, the user can always provide a base ring
by passing it in as the first argument to the `Singular.jl` function. Note that
if the first argument to the `Singular.jl` version of a library procedure is a
polynomial or non-commutative ring, then this is automatically assumed to be the
base ring. Hence, if a procedure in the Singular language takes a ring as a
first argument, you will have to pass that ring as the second argument after
specifying the base ring in the first argument.

**Examples**

This example illustrates passing Singular lists and providing the base ring.

```jldoctest
julia> r0, (x, y, z, t) = polynomial_ring(QQ, ["x", "y", "z", "t"], ordering=ordering_lp());

julia> Singular.LibGeneral.sort([x, y])
ERROR: `intvec` may be passed in as Vector{Int}. All other vectors (`list` in Singular) must be passed in as Vector{Any} along with an explicit base ring in the first argument

julia> Singular.LibGeneral.sort(r0, Any[x, y])
2-element Vector{Vector}:
 spoly{n_Q}[y, x]
 [2, 1]
```

This example illustrates the base ring inference:

```jldoctest
julia> AA, (x, y, z, t) = polynomial_ring(QQ, ["x", "y", "z", "t"]);

julia> D = zero_matrix(AA, 4, 4);

julia> D[1,2] = -z; D[1,3] = 2*x; D[2,3] = -2*y;

julia> A, (x, y, z, t) = GAlgebra(AA, 1, D);

julia> Singular.LibNctools.isCentral(x)   # base ring A is inferred from x
0

julia> Singular.LibCentral.center(A, 3)   # base ring cannot be inferred from the plain Int 3
Singular ideal over Singular G-Algebra (QQ),(x,y,z,t),(dp(4),C) with generators (t, 4*x*y + z^2 - 2*z)
```

## Global Interpreter Variables

The function `Singular.call_interpreter` can be used to execute arbitrary
strings inside the Singular interpreter, and the function
`Singular.lookup_library_symbol` fetches results from the interpreter.

```@docs
Singular.lookup_library_symbol(package::String, name::String)
```

## Global Kernel Variables

The global variables `degBound` and `multBound` can be used in a local fashion.
As with any global variable, their usage should be accompanied with caution.

```@docs
with_degBound(f, degb::Integer)
```

```@docs
with_multBound(f, degb::Integer)
```

The following [options](https://www.singular.uni-kl.de/Manual/4-3-0/sing_318.htm#SEC358)
are available. The usage of, say, the option `infRefTail`
would be as `with_infRefTail(f, flag::Bool)` where the same do-block syntax
can be used as with the degree bounds.

```
fastHC, infRedTail, lazy, length, notBuckets, prot, qringNF, redTail, redThrough
```

**Examples**

```jldoctest
julia> r, (x,y,z) = polynomial_ring(QQ, ["x", "y", "z"], ordering=ordering_ds());

julia> i = Ideal(r, [x^7+y^7+z^6,x^6+y^8+z^7,x^7+y^5+z^8,x^2*y^3+y^2*z^3+x^3*z^2,x^3*y^2+y^3*z^2+x^2*z^3]);

julia> degree(std(i))   # default behaviour of no multiplicity bound
(0, 86)

julia> with_multBound(100) do
           # run with a multiplicity bound of 100
           return degree(std(i))
       end
(0, 98)

julia> degree(std(i))   # back to default behaviour
(0, 86)

julia> gens(std(i))
11-element Vector{spoly{n_Q}}:
 x^3*y^2 + y^3*z^2 + x^2*z^3
 x^2*y^3 + x^3*z^2 + y^2*z^3
 y^5 + x^7 + z^8
 x^6 + z^7 + y^8
 x^4*z^2 - y^4*z^2 - x^2*y*z^3 + x*y^2*z^3
 z^6 + x^7 + y^7
 y^4*z^3 - y^3*z^4 - x^2*z^5 - x^9
 x^3*y*z^4 - x^2*y^2*z^4 + x*y^3*z^4 - y^4*z^4 + x^3*z^5 - x^2*y*z^5
 x^3*z^5
 x^2*y*z^5 + y^3*z^5 + x^2*z^6
 x*y^3*z^5

julia> gens(with_degBound(5) do; return std(i); end)
5-element Vector{spoly{n_Q}}:
 x^3*y^2 + y^3*z^2 + x^2*z^3
 x^2*y^3 + x^3*z^2 + y^2*z^3
 y^5 + x^7 + z^8
 x^6 + z^7 + y^8
 z^6 + x^7 + y^7

julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> with_prot(true) do; return std(Ideal(R, x^5 - y*x + 1, y^6*x + x^2 + y^3)); end
[4294967295:2]5s7s11s1214-s15
product criterion:1 chain criterion:1
Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x^5 - x*y + 1, x*y^6 + y^3 + x^2, x^4*y^3 - y^6 - y^4 - x, y^9 + y^7 + x^3*y^3 + x*y^3 + x*y - 1)
```

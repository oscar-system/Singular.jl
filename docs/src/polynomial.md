```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Multivariate polynomials

Singular.jl allows the creation of multivariate polynomials over any of the coefficient
rings described above.

The default multivariate polynomial type in Singular.jl is the Singular `spoly` type.

The associated polynomial ring is represented by a parent object which can be
constructed by a call to the `polynomial_ring` constructor.

The types of the polynomial ring parent objects and elements thereof are given in the
following table according to the library providing them.

 Library        | Element type  | Parent type
----------------|---------------|-----------------------
Singular        | `spoly{T}`    | `Singular.PolyRing{T}`

These types are parameterised by the type of elements in the coefficient ring of the
polynomials.

All polynomial types belong directly to the abstract type `MPolyRingElem` and
all the polynomial ring parent object types belong to the abstract type `MPolyRing`.

## Multivariate polynomial functionality

Singular.jl polynomials provides all the Multivariate Polynomial Ring functionality
described by AbstractAlgebra.jl.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/mpolynomial>

In particular, Singular polynomials are sparse distributed, but do not have random
access. Instead, they implement iterator access to terms. This is due to their storage
in a linked list, for efficient implementation of Groebner basis algorithms.

Some polynomial rings may also implement part of the Euclidean Ring interface, where
this is appropriate.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/euclidean_interface>

Below, we describe the functionality that is specific to the Singular multivariate
polynomials that is not documented in the general multivariate interface.

### Constructors

```julia
polynomial_ring(R::Union{Ring, Field}, s::AbstractVector{<:VarName};
               cached::Bool = true, ordering = :degrevlex,
               ordering2::Symbol = :comp1min, degree_bound::Int = 0)
```

Returns a tuple, $S, x$ consisting of a multivariate polynomial ring $S$ and an array $x$
of variables (from which polynomials can be constructed). The ring $R$ must be a valid
Singular coefficient ring, or any Nemo/AbstractAlgebra coefficient ring. The array $s$
must be a list of strings corresponding to how the variables will be printed. By default,
there will only be one Singular polynomial ring in the system for each combination of
coefficient ring, list of variable names, ordering and degree bound. This is
accomplished by making use of a global cache. If this is not the desired behaviour,
`false` can be passed to the optional argument `cached`.

If the first ordering `ordering` is specified as a symbol, then
two orderings can be specified, one for term ordering of the polynomials, and another
for ordering of module components. They can occur in either order, the first taking
precedence over the other, when the polynomials are used to represent module generators.
If either is not specified, the indicated default is used.
The options for polynomial term ordering are, `:lex`, `:deglex`, `:degrevlex`,
`:neglex`, `:negdeglex` and `:negdegrevlex`, and the options for module component ordering
are `comp1min` and `comp1max`.

If the first ordering `ordering` is specified as a non-symbol, the second ordering
`ordering2` will be ignored. For specifying non-symbolic term orderings, please
see the Term orderings section below.

If one has an a priori bound on the degree in each variable of a polynomial (including
for all intermediate computations in this ring), one can specify it using the
`degree_bound` optional parameter. Singular may then be able to use a more efficient
representation internally, which will use less memory and allow for faster arithmetic.
By default, Singular uses a bound of 16 bits internally for the exponent of each
variable, however this is a signed value, so that the default is for nonnegative
exponents that fit in 15 bits.

Note that internally, Singular may use a higher bound than specified, if it will not
increase the amount of storage required.

**Examples**

```jldoctest
julia> R, (x, y, z) = polynomial_ring(ZZ, ["x", "y", "z"])
(Singular polynomial ring (ZZ),(x,y,z),(dp(3),C), spoly{n_Z}[x, y, z])

julia> S, vars = polynomial_ring(QQ, ["x", "y"]; ordering=:deglex)
(Singular polynomial ring (QQ),(x,y),(Dp(2),C), spoly{n_Q}[x, y])

julia> T, x = polynomial_ring(ZZ, ["x$i" for i in 1:5];
              ordering=:comp1max, ordering2=:degrevlex, degree_bound=5)
(Singular polynomial ring (ZZ),(x1,x2,x3,x4,x5),(c,dp(5),L(5)), spoly{n_Z}[x1, x2, x3, x4, x5])
```

See also the convenience macros below for simple use cases.

The following function allows creating a Singular polynomial ring from a given
polynomial ring of type AbstractAlgebra.Generic.MPolyRing:

```julia
polynomial_ring(R::AbstractAlgebra.Generic.MPolyRing{T}; cached::Bool = true, ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min, degree_bound::Int = 0)  where {T <: RingElement}
```

Polynomials can be constructed using arithmetic operations on the generators,
but this can be extremely inefficient. For this purpose, Singular polynomials
support construction using a build context, as described in the AbstractAlgebra
documentation:

<https://nemocas.github.io/AbstractAlgebra.jl/latest/mpolynomial>

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Singular polynomial ring (ZZ),(x,y),(dp(2),C), spoly{n_Z}[x, y])

julia> C = MPolyBuildCtx(R)
Builder for an element of Singular polynomial ring (ZZ),(x,y),(dp(2),C)

julia> push_term!(C, ZZ(1), [1, 2])
x*y^2

julia> push_term!(C, ZZ(3), [1, 1])
x*y^2 + 3*x*y

julia> push_term!(C, -ZZ(1), [0, 1])
x*y^2 + 3*x*y - y

julia> f = finish(C)
x*y^2 + 3*x*y - y
```

### Term orderings

A general term ordering can be constructed as a product of one or more of the
following block orderings.

```@docs
ordering_lp(nvars::Int = 1)
```

```@docs
ordering_rp(nvars::Int = 1)
```

```@docs
ordering_ip(nvars::Int = 1)
```

```@docs
ordering_dp(nvars::Int = 1)
```

```@docs
ordering_Dp(nvars::Int = 1)
```

```@docs
ordering_wp(w::Vector{Int})
```

```@docs
ordering_Wp(w::Vector{Int})
```

```@docs
ordering_Ip(nvars::Int = 1)
```

```@docs
ordering_ls(nvars::Int = 1)
```

```@docs
ordering_rs(nvars::Int = 1)
```

```@docs
ordering_is(nvars::Int = 1)
```

```@docs
ordering_ds(nvars::Int = 1)
```

```@docs
ordering_Ds(nvars::Int = 1)
```

```@docs
ordering_ws(w::Vector{Int})
```

```@docs
ordering_Ws(w::Vector{Int})
```

```@docs
ordering_a(w::Vector{Int})
```

```@docs
ordering_M(m::Matrix{Int}; checked::Bool = true)
```

```@docs
ordering_C()
```

```@docs
ordering_c()
```

**Examples**

```julia
polynomial_ring(QQ, "x".*string.(1:8), ordering = ordering_M([1 2; 3 5])*ordering_lp(3)*ordering_wp([1, 2, 3]))

polynomial_ring(QQ, "x".*string.(1:5), ordering = ordering_dp(3)*ordering_dp())
```

### Polynomial ring macros

For convenience, we provide some macros for constructing polynomial rings and injecting
the variables into scope. These are easier to use, but have some limitations. In
particular, they can only be used at the top level by the user, and cannot be used
programmatically or in library code (it is not possible to inject an arbitrary number
of variables into scope inside a function).

The macros are designed for simple use cases, and do not offer the full power of the
most general constructor above.

```julia
@polynomial_ring(R, s, n, o)
```

Given a coefficient ring $R$, a root variable name, e.g. `"x"`, a number of variable
$n$ and a polynomial term ordering `o`, create the variables `x1, x2, ..., xn` and
inject them into scope, and return the corresponding polynomial ring `S`.

```julia
@polynomial_ring(R, s, n)
```

As per the previous macro, with a default of `:degrevlex` for the polynomial term
ordering.

**Examples**

```jldoctest
julia> S = @polynomial_ring(ZZ, "x", 5, :deglex)
Singular polynomial ring (ZZ),(x1,x2,x3,x4,x5),(Dp(5),C)

julia> T = @polynomial_ring(QQ, "y", 10)
Singular polynomial ring (QQ),(y1,y2,y3,y4,y5,y6,y7,y8,y9,y10),(dp(10),C)
```

### Basic manipulation

```@docs
has_global_ordering(R::PolyRing)
```

```@docs
has_mixed_ordering(R::PolyRing)
```

```@docs
has_local_ordering(R::PolyRing)
```

```@docs
degree_bound(R::PolyRing)
```

```@docs
total_degree(p::spoly)
```

```@docs
order(p::spoly)
```

**Examples**

```jldoctest
julia> R = @polynomial_ring(ZZ, "x", 3)
Singular polynomial ring (ZZ),(x1,x2,x3),(dp(3),C)

julia> n = nvars(R)
3

julia> has_global_ordering(R) == true
true

julia> c = characteristic(R)
0

julia> L = degree_bound(R)
1048575

julia> exps = leading_exponent_vector(x1*x2 + 3x1*x2^2 + x3 + 2)
3-element Vector{Int64}:
 1
 2
 0

julia> deg = total_degree(x1*x2 + 3x1*x2^2 + x3 + 2)
3

julia> ord = order(x1*x2 + 3x1*x2^2 + x3 + 2)
0
```

### Differential functions

Working over any coefficient ring, basic functionality involving
differential operations is available.

```@docs
jet(::spoly{T}, ::Int) where T <: Nemo.RingElem
```

```@docs
derivative(::spoly{T}, ::Int) where T <: Nemo.RingElem
```

```@docs
derivative(::spoly{T}, ::spoly{T}) where T <: Nemo.RingElem
```

```@docs
jacobian_ideal(::spoly{T}) where T <: Nemo.RingElem
```

```@docs
jacobian_matrix(p::spoly{T}) where T <: Nemo.RingElem
```

```@docs
jacobian_matrix(A::Vector{spoly{T}}) where T <: Nemo.RingElem
```

**Examples**

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Singular polynomial ring (QQ),(x,y,z),(dp(3),C), spoly{n_Q}[x, y, z])

julia> f = x^2*y*z + z^2*x + x*y*z
x^2*y*z + x*y*z + x*z^2

julia> g = jet(f, 3)
x*y*z + x*z^2

julia> derivative(f, 1)
2*x*y*z + y*z + z^2

julia> derivative(f, y)
x^2*z + x*z

julia> J = jacobian_ideal(f)
Singular ideal over Singular polynomial ring (QQ),(x,y,z),(dp(3),C) with generators (2*x*y*z + y*z + z^2, x^2*z + x*z, x^2*y + x*y + 2*x*z)

julia> Jf1 = jacobian_matrix(f)
[2*x*y*z + y*z + z^2
x^2*z + x*z
x^2*y + x*y + 2*x*z]

julia> Jf2 = jacobian_matrix([f, g])
[2*x*y*z + y*z + z^2, x^2*z + x*z, x^2*y + x*y + 2*x*z
y*z + z^2, x*z, x*y + 2*x*z]
```

### Content and primitive part

When coefficient rings have a meaningful GCD function, the following functions are
available.

```@docs
Singular.primpart(x::spoly)
```

```@docs
Singular.content(x::spoly)
```

**Examples**

```jldoctest
julia> R = @polynomial_ring(ZZ, "x", 2)
Singular polynomial ring (ZZ),(x1,x2),(dp(2),C)

julia> f = 3x1^2 + 3x1*x2 + 6x2^2
3*x1^2 + 3*x1*x2 + 6*x2^2

julia> p = primpart(f)
x1^2 + x1*x2 + 2*x2^2

julia> c = content(f)
3
```

### Homogeneous polynomials

```@docs
homogenize(p::spoly{T}, v::spoly{T}) where T <: Nemo.RingElem
```

### Multivariate Factorisation

For the Singular base fields `QQ` and `Fp` a function to compute a
squarefree factorization is available.

**Examples**

```julia
julia> R = @polynomial_ring(QQ, "x", 4)
Singular Polynomial Ring (QQ),(x1,x2,x3,x4),(dp(4),C)

julia> f = 123*(57*x2^3 + x4^5)^3*(x1^2 + x1+1)^2*(x1 + x2*x3)^2;

julia> Fac = factor(f)
123 * (x4^5 + 57*x2^3)^3 * (x1^2 + x1 + 1)^2 * (x2*x3 + x1)^2
```

For the Singular base rings `QQ`, `ZZ` and `Fp` a function to compute the
multivariate factorization is available.

**Examples**

```julia
julia> R = @polynomial_ring(ZZ, "x", 4)
Singular Polynomial Ring (ZZ),(x1,x2,x3,x4),(dp(4),C)

julia> f = 123*(57*x2^3 + x4^5)^3*(x1^2 + x1+1)^2*(x1 + x2*x3)^2;

julia> Fac = factor(f)
123 * (x2*x3 + x1)^2 * (x4^5 + 57*x2^3)^3 * (x1^2 + x1 + 1)^2
```

### Change of coefficient rings

It is possible to change the coefficient ring of a given polynomial $p$ via
the function 'change_base_ring'.

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Singular polynomial ring (ZZ),(x,y),(dp(2),C), spoly{n_Z}[x, y])

julia> p = x^5 + y^3+1
x^5 + y^3 + 1

julia> p2 = change_base_ring(QQ, p)
x^5 + y^3 + 1

julia> parent(p2)
Singular polynomial ring (QQ),(x,y),(dp(2),C)
```

It also possible to work with Nemo rings by casting to a suitable Singular type
via `CoefficientRing`.

**Examples**

```julia
R, (x, y) = polynomial_ring(ZZ, ["x", "y"])

p = x^5 + y^3+1

p2 change_base_ring(CoefficientRing(Nemo.QQ), p)
```

```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Ideals

Singular.jl allows the creation of ideals over a Singular polynomial ring. These are
internally stored as a list of (polynomial) generators. This list of generators can also
have the property of being a Groebner basis.

The default ideal type in Singular.jl is the Singular `sideal` type.

Ideals objects have a parent object which represents the set of ideals they belong to,
the data for which is given by the polynomial ring their generators belong to.

The types of ideals and associated parent objects are given in the following table
according to the library providing them.

 Library        | Element type  | Parent type
----------------|---------------|-----------------------
Singular        | `sideal{T}`    | `Singular.IdealSet{T}`

These types are parameterised by the type of elements of the polynomial ring over which
the ideals are defined.

All ideal types belong directly to the abstract type `Module{T}` and
all the ideal set parent object types belong to the abstract type `Set`.

## Ideal functionality

Singular.jl ideals implement standard operations one would expect on modules.
These include:

 * Operations common to all AbstractAlgebra objects, such as `parent`, `base_ring`,
   `elem_type`, `parent_type`, `parent`, `deepcopy`, etc.

 * Addition

Also implements is the following operations one expects for ideals:

 * Multiplication

 * Powering

Below, we describe all of the functionality for Singular.jl ideals that is not included
in this list of basic operations.

### Constructors

Given a Singular polynomial ring $R$, the following constructors are available for
creating ideals.

```julia
Ideal(R::PolyRing{T}, ids::spoly{T}...) where T <: Nemo.RingElem
Ideal(R::PolyRing{T}, ids::Vector{spoly{T}}) where T <: Nemo.RingElem
```

Construct the ideal over the polynomial ring $R$ whose (polynomial) generators are given
by the given parameter list or array of polynomials, respectively. The list may be
empty, resulting in the zero ideal.

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Singular Polynomial Ring (ZZ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Z}[x, y])

julia> I1 = Ideal(R, x*y + 1, x^2)
Singular ideal over Singular Polynomial Ring (ZZ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x*y + 1, x^2)

julia> I2 = Ideal(R, [x*y + 1, x^2])
Singular ideal over Singular Polynomial Ring (ZZ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x*y + 1, x^2)
```

### Basic manipulation

```@docs
ngens(::sideal)
```

```@docs
gens(::sideal)
```

Singular.jl overloads the `setindex!` and `getindex` functions so that one can access
the generators of an ideal using array notation.

```julia
I[n::Int]
```

```@docs
iszero(::sideal)
```

```@docs
is_zerodim(I::sideal)
```

```@docs
dimension(I::sideal{spoly{T}}) where T <: Nemo.RingElem
```

```@docs
is_constant(::sideal)
```

```@docs
is_var_generated(::sideal)
```

```@docs
normalize!(::sideal)
```

```@docs
interreduce(I::sideal{S}) where {T <: Nemo.RingElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
(Singular Polynomial Ring (ZZ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Z}[x, y])

julia> I = Ideal(R, x^2 + 1, x*y)
Singular ideal over Singular Polynomial Ring (ZZ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + 1, x*y)

julia> n = ngens(I)
2

julia> p = I[1]
x^2 + 1

julia> I[1] = 2x + y^2
y^2 + 2*x

julia> is_constant(I) == false
true

julia> is_var_generated(I) == false
true

julia> is_zerodim(I) == false
ERROR: Not a Groebner basis

julia> S, (u, v) = polynomial_ring(QQ, ["u", "v"])
(Singular Polynomial Ring (QQ),(@OSCAR@u@1,@OSCAR@v@1),(dp(2),C), spoly{n_Q}[u, v])

julia> J = Ideal(S, u^2 + 1, u*v)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@u@1,@OSCAR@v@1),(dp(2),C) with generators (u^2 + 1, u*v)

julia> dimension(std(J)) == 0
true
```

### Containment

```@docs
contains{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})
```

**Examples**

```jldoctest
julia> R, (x , y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2 + 1, x*y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + 1, x*y)

julia> J = Ideal(R, x^2 + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + 1)

julia> contains(I, J) == true
true
```

### Comparison

Checking whether two ideals are algebraically equal is very expensive, as it usually
requires computing Groebner bases. Therefore we do not overload the `==` operator for
ideals. Instead we have the following two functions.

```@docs
isequal(I1::sideal{S}, I2::sideal{S}) where S <: SPolyUnion
```

```@docs
equal(I1::sideal{S}, I2::sideal{S}) where S <: SPolyUnion
```

**Examples**

```jldoctest
julia> R, (x , y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2 + 1, x*y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + 1, x*y)

julia> J = Ideal(R, x^2 + x*y + 1, x^2 - x*y + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + x*y + 1, x^2 - x*y + 1)

julia> isequal(I, J) == false
true

julia> equal(I, J) == true
true
```

### Intersection

```@docs
intersection(I::sideal{S}, J::sideal{S}) where {T <: Nemo.RingElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```jldoctest
julia> R, (x , y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2 + 1, x*y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + 1, x*y)

julia> J = Ideal(R, x^2 + x*y + 1, x^2 - x*y + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + x*y + 1, x^2 - x*y + 1)

julia> V = intersection(I, J)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (y, x^2 - x*y + 1)
```

### Quotient

```@docs
quotient(I::sideal{S}, J::sideal{S}) where S <: spoly
```

```@docs
quotient(I::sideal{S}, J::sideal{S}) where S <: spluralg
```

**Examples**

```jldoctest
julia> R, (x , y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2 + 1, x*y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + 1, x*y)

julia> J = Ideal(R, x + y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x + y)

julia> V = quotient(I, J)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (y, x^2 + 1)
```

### Leading terms

```@docs
lead(I::sideal{S}) where S <: SPolyUnion
```

**Examples**

```jldoctest
julia> R, (x , y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2 + 1, x*y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + 1, x*y)

julia> V = lead(I)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2, x*y)
```

### Homogeneous ideals

```@docs
is_homogeneous(I::sideal)
```

```@docs
homogenize(I::sideal{S}, v::S) where S <: spoly
```

### Saturation

```@docs
saturation(I::sideal{T}, J::sideal{T}) where T <: Nemo.RingElem
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, (x^2 + x*y + 1)*(2y^2+1)^3, (2y^2 + 3)*(2y^2+1)^2)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (8*x^2*y^6 + 8*x*y^7 + 12*x^2*y^4 + 12*x*y^5 + 8*y^6 + 6*x^2*y^2 + 6*x*y^3 + 12*y^4 + x^2 + x*y + 6*y^2 + 1, 8*y^6 + 20*y^4 + 14*y^2 + 3)

julia> J = Ideal(R, 2y^2 + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (2*y^2 + 1)

julia> S = saturation(I, J)
(Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (2*y^2 + 3, x^2 + x*y + 1), 2)
```

### Standard basis

```@docs
std(I::sideal{S}; complete_reduction::Bool=false) where S <: SPolyUnion
```

```@docs
fglm(I::sideal{spoly{T}}, ordering::Symbol) where T <: Nemo.RingElem
```

```@docs
satstd(I::sideal{spoly{T}}, J::sideal{spoly{T}}) where T <: Nemo.FieldElem
```

```@docs
lift_std(M::sideal{S}; complete_reduction::Bool = false) where S <: spoly
```

```@docs
lift_std_syz(M::sideal{S}; complete_reduction::Bool = false) where S <: spoly
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2 + x*y + 1, 2*y^2 + 3)

julia> J = Ideal(R, 2*y^2 + 3, x^2 + x*y + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (2*y^2 + 3, x^2 + x*y + 1)

julia> A = std(I)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (2*y^2 + 3, x^2 + x*y + 1)

julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, (x*y + 1)*(2x^2*y^2 + x*y - 2) + 2x*y^2 + x, 2x*y + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (2*x^3*y^3 + 3*x^2*y^2 + 2*x*y^2 - x*y + x - 2, 2*x*y + 1)

julia> J = Ideal(R, x)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x)

julia> B = satstd(I, J)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x - y - 1, 2*y^2 + 2*y + 1)

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering = :lex)
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@z@1),(lp(3),C), spoly{n_Q}[x, y, z])

julia> I = Ideal(R, y^3+x^2, x^2*y+x^2, x^3-x^2, z^4-x^2-y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@z@1),(lp(3),C) with generators (x^2 + y^3, x^2*y + x^2, x^3 - x^2, -x^2 - y + z^4)

julia> J = fglm(I, :degrevlex)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@z@1),(lp(3),C) with generators (z^12, y*z^4 - z^8, y^2 + y - z^8 - z^4, x*y - x*z^4 - y + z^4, x^2 + y - z^4)
```

### Reduction

```@docs
reduce(I::sideal{S}, G::sideal{S}) where S <: SPolyUnion
```

```@docs
reduce(p::S, G::sideal{S}) where S <: SPolyUnion
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> f = x^2*y + 2y + 1
x^2*y + 2*y + 1

julia> g = y^2 + 1
y^2 + 1

julia> I = Ideal(R, (x^2 + 1)*f + (x + y)*g + x + 1, (2y^2 + x)*f + y)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^4*y + 3*x^2*y + x*y^2 + y^3 + x^2 + 2*x + 3*y + 2, 2*x^2*y^3 + x^3*y + 4*y^3 + 2*x*y + 2*y^2 + x + y)

julia> J = std(Ideal(R, f, g))
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (y^2 + 1, x^2 - y + 2)

julia> V = reduce(I, J)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x + 1, y)

julia> h1 = (x^2 + 1)*f + (x + y)*g + x + 1
x^4*y + 3*x^2*y + x*y^2 + y^3 + x^2 + 2*x + 3*y + 2

julia> h2 = reduce(h1, J)
x + 1
```

### Elimination

```@docs
eliminate(I::sideal{S}, polys::S...) where {T <: Nemo.RingElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```jldoctest
julia> R, (x, y, t) = polynomial_ring(QQ, ["x", "y", "t"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@t@1),(dp(3),C), spoly{n_Q}[x, y, t])

julia> I = Ideal(R, x - t^2, y - t^3)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@t@1),(dp(3),C) with generators (-t^2 + x, -t^3 + y)

julia> J = eliminate(I, t)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@t@1),(dp(3),C) with generators (x^3 - y^2)
```

### Syzygies

```@docs
syz(::sideal)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2*y + 2*y + 1, y^2 + 1)

julia> F = syz(I)
Singular Module over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), with Generators:
@OSCAR@x@1^2*@OSCAR@y@1*gen(2)-@OSCAR@y@1^2*gen(1)+2*@OSCAR@y@1*gen(2)+gen(2)-gen(1)

julia> M = Singular.Matrix(I)
[x^2*y + 2*y + 1, y^2 + 1]

julia> N = Singular.Matrix(F)
[-y^2 - 1
x^2*y + 2*y + 1]

julia> iszero(M*N)  # check they are actually syzygies
true
```

### Free resolutions

```@docs
fres{T <: Nemo.FieldElem}(::sideal{spoly{T}}, ::Int, ::String)
```

```@docs
sres{T <: Nemo.FieldElem}(::sideal{spoly{T}}, ::Int)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(dp(2),C) with generators (x^2*y + 2*y + 1, y^2 + 1)

julia> F1 = fres(std(I), 0)
Singular Resolution:
R^1 <- R^2 <- R^1

julia> F2 = sres(std(I), 2)
Singular Resolution:
R^1 <- R^2 <- R^1
```

### Differential operations

```@docs
jet(I::sideal{S}, n::Int) where {T <: Nemo.RingElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@z@1),(dp(3),C), spoly{n_Q}[x, y, z])

julia> I = Ideal(R, x^5 - y^2, y^3 - x^6 + z^3)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@z@1),(dp(3),C) with generators (x^5 - y^2, -x^6 + y^3 + z^3)

julia> J1 = jet(I, 3)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@z@1),(dp(3),C) with generators (-y^2, y^3 + z^3)
```

### Operations on zero-dimensional ideals

```@docs
vdim(I::sideal{S}) where {T <: Nemo.FieldElem, S <: Union{spoly{T}, spluralg{T}}}
```

```@docs
kbase(I::sideal{S}) where {T <: Nemo.FieldElem, S <: Union{spoly{T}, spluralg{T}}}
```

```@docs
highcorner(I::sideal{S}) where {T <: Nemo.FieldElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]; ordering=:negdegrevlex)
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(ds(2),C), spoly{n_Q}[x, y])

julia> I = Ideal(R, 3*x^2 + y^3, x*y^2)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(ds(2),C) with generators (3*x^2 + y^3, x*y^2)

julia> I = std(I)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(ds(2),C) with generators (3*x^2 + y^3, x*y^2, y^5)

julia> n = vdim(I)
7

julia> J = kbase(I)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(ds(2),C) with generators (y^4, y^3, y^2, x*y, y, x, 1)

julia> f = highcorner(I)
y^4
```

### Operations over local rings

If the base ring `R` is a local ring, a minimal generating set can be computed
using the following function

```@docs
minimal_generating_set(I::sideal{S}) where S <: spoly
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]; ordering=:negdegrevlex)
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(ds(2),C), spoly{n_Q}[x, y])

julia> has_local_ordering(R) == true
true

julia> I = Ideal(R, y, x^2, (1 + y^3) * (x^2 - y))
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1),(ds(2),C) with generators (y, x^2, -y + x^2 - y^4 + x^2*y^3)

julia> min = minimal_generating_set(I)
2-element Vector{spoly{n_Q}}:
 x^2
 y
```

### Independent sets of monomial ideals

Let $I$ be an ideal of $K[x_1, ..., x_n].$ An `independent set` is a
subset $u \subseteq \{x_1, ..., x_n\},$ such that $I \cap K[u]= 0.$ In case
$u$ cannot be enlarged, it is called `non-extendable independent set`.
If in addition $|u| = dim(K[x_1, ..., x_n]/I),$ $u$ is called
`maximal independent set`.
Using Singular.jl one can compute non-extendable, resp. maximal independent
sets for monomial ideals.
If an arbitrary ideal $I$ is passed to the function, the computation is performed on
the leading ideal of $I$.

```@docs
independent_sets(I::sideal{spoly{T}}) where T <: Nemo.FieldElem
```

```@docs
maximal_independent_set(I::sideal{spoly{T}}; all::Bool = false) where T <: Nemo.FieldElem
```

```jldoctest
julia> R, (x, y, u, v, w) = polynomial_ring(QQ, ["x", "y", "u", "v", "w"])
(Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@u@1,@OSCAR@v@1,@OSCAR@w@1),(dp(5),C), spoly{n_Q}[x, y, u, v, w])

julia> has_local_ordering(R) == true
false

julia> I = Ideal(R, x*y*w, y*v*w, u*y*w, x*v)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@u@1,@OSCAR@v@1,@OSCAR@w@1),(dp(5),C) with generators (x*y*w, y*v*w, y*u*w, x*v)

julia> I = std(I)
Singular ideal over Singular Polynomial Ring (QQ),(@OSCAR@x@1,@OSCAR@y@1,@OSCAR@u@1,@OSCAR@v@1,@OSCAR@w@1),(dp(5),C) with generators (x*v, y*v*w, y*u*w, x*y*w)

julia> L1 = independent_sets(I)
5-element Vector{Vector{spoly{n_Q}}}:
 [x, y, u]
 [y, u, v]
 [x, u, w]
 [u, v, w]
 [y, w]

julia> L2 = maximal_independent_set(I)
3-element Vector{spoly{n_Q}}:
 x
 y
 u

julia> L3 = maximal_independent_set(I, all = true)
4-element Vector{Vector{spoly{n_Q}}}:
 [x, y, u]
 [y, u, v]
 [x, u, w]
 [u, v, w]
```

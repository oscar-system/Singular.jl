```@meta
CurrentModule = Singular
```

# Ideals

Singular.jl allows the creation of ideals over a Singular polynomial ring. These are
internally stored as a list of (polynomial) generators. This list of generators can also
have the property of being a Groebner basis.

The default ideal type in Singular.jl is the Singular `sideal` type.

Ideals objects have a parent object which represents the set of ideals they belong to,
the data for which is given by the polynomial ring their generators belong to.

The types of ideals and associated parent objects are given in the following table
according to the library provding them.

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

```julia
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

I1 = Ideal(R, x*y + 1, x^2)
I2 = Ideal(R, [x*y + 1, x^2])
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

```
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)

n = ngens(I)
p = I[1]
I[1] = 2x + y^2
is_constant(I) == false
is_var_generated(I) == false
is_zerodim(I) == false

S, (u, v) = PolynomialRing(QQ, ["u", "v"])
J = Ideal(S, u^2 + 1, u*v)
dimension(std(J)) == 0
```

### Containment

```@docs
contains{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})
```

**Examples**

```julia
R, (x , y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)
J = Ideal(R, x^2 + 1)

contains(I, J) == true
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

```julia
R, (x , y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)
J = Ideal(R, x^2 + x*y + 1, x^2 - x*y + 1)

isequal(I, J) == false
equal(I, J) == true
```

### Intersection

```@docs
intersection(I::sideal{S}, J::sideal{S}) where {T <: Nemo.RingElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```julia
R, (x , y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)
J = Ideal(R, x^2 + x*y + 1, x^2 - x*y + 1)

V = intersection(I, J)
```

### Quotient

```@docs
quotient(I::sideal{S}, J::sideal{S}) where S <: spoly
```

```@docs
quotient(I::sideal{S}, J::sideal{S}) where S <: spluralg
```

**Examples**

```julia
R, (x , y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)
J = Ideal(R, x + y)

V = quotient(I, J)
```

### Leading terms

```@docs
lead(I::sideal{S}) where S <: SPolyUnion
```

**Examples**

```julia
R, (x , y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)

V = lead(I)
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

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, (x^2 + x*y + 1)*(2y^2+1)^3, (2y^2 + 3)*(2y^2+1)^2)
J = Ideal(R, 2y^2 + 1)

S = saturation(I, J)
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

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
J = Ideal(R, 2*y^2 + 3, x^2 + x*y + 1)

A = std(I)

R, (x, y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, (x*y + 1)*(2x^2*y^2 + x*y - 2) + 2x*y^2 + x, 2x*y + 1)
J = Ideal(R, x)

B = satstd(I, J)

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering = :lex)
I = Ideal(R, y^3+x^2, x^2*y+x^2, x^3-x^2, z^4-x^2-y)
J = fglm(I, :degrevlex)
```

### Reduction

```@docs
reduce(I::sideal{S}, G::sideal{S}) where S <: SPolyUnion
```

```@docs
reduce(p::S, G::sideal{S}) where S <: SPolyUnion
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

f = x^2*y + 2y + 1
g = y^2 + 1

I = Ideal(R, (x^2 + 1)*f + (x + y)*g + x + 1, (2y^2 + x)*f + y)
J = std(Ideal(R, f, g))

V = reduce(I, J)

h1 = (x^2 + 1)*f + (x + y)*g + x + 1

h2 = reduce(h1, J)
```

### Elimination

```@docs
eliminate(I::sideal{S}, polys::S...) where {T <: Nemo.RingElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```julia
R, (x, y, t) = PolynomialRing(QQ, ["x", "y", "t"])

I = Ideal(R, x - t^2, y - t^3)

J = eliminate(I, t)
```

### Syzygies

```@docs
syz(::sideal)
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)

F = syz(I)

M = Singular.Matrix(I)
N = Singular.Matrix(F)

# check they are actually syzygies
iszero(M*N)
```

### Free resolutions

```@docs
fres{T <: Nemo.FieldElem}(::sideal{spoly{T}}, ::Int, ::String)
```

```@docs
sres{T <: Nemo.FieldElem}(::sideal{spoly{T}}, ::Int)
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)

F1 = fres(std(I), 0)
F2 = sres(std(I), 2)
```

### Differential operations

```@docs
jet(I::sideal{S}, n::Int) where {T <: Nemo.RingElem, S <: Union{spoly{T}, spluralg{T}}}
```

**Examples**

```julia
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

I = Ideal(R, x^5 - y^2, y^3 - x^6 + z^3)

J1 = jet(I, 3)
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

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"]; ordering=:negdegrevlex)

I = Ideal(R, 3*x^2 + y^3, x*y^2)

I = std(I)

n = vdim(I)
J = kbase(I)
f = highcorner(I)
```

### Operations over local rings

If the base ring `R` is a local ring, a minimal generating set can be computed
using the following function

```@docs
minimal_generating_set(I::sideal{S}) where S <: spoly
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"]; ordering=:negdegrevlex)

has_local_ordering(R) == true

I = Ideal(R, y, x^2, (1 + y^3) * (x^2 - y))

min = minimal_generating_set(I)
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

```julia
R, (x, y, u, v, w) = PolynomialRing(QQ, ["x", "y", "u", "v", "w"])

has_local_ordering(R) == true

I = Ideal(R, x*y*w, y*v*w, u*y*w, x*v)

I = std(I)

L1 = independent_sets(I)

L2 = maximal_independent_set(I)

L3 = maximal_independent_set(I, all = true)

```

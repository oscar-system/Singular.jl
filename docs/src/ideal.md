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
Ideal(R::PolyRing{T}, ids::Array{spoly{T}, 1}) where T <: Nemo.RingElem
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

Singular.jl overloads the `setindex!` and `getindex` functions so that one can access
the generators of an ideal using array notation.

```julia
I[n::Int]
```

```@docs
iszero(::sideal)
```

```@docs
iszerodim(::sideal)
```

```@docs
isconstant(::sideal)
```

```@docs
isvar_generated(::sideal)
```

```@docs
normalize!(::sideal)
```

```@docs
radical(::sideal)
```

**Examples**

```
R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)

n = ngens(I)
p = I[1]
I[1] = 2x + y^2
isconstant(I) == false
isvar_generated(I) == false
iszerodim(I) == false
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
isequal{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})
```

```@docs
equal{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})
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
intersection{T <: Nemo.RingElem}(::sideal{T}, ::sideal{T})
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
quotient{T <: Nemo.RingElem}(::sideal{T}, ::sideal{T})
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
lead(::sideal)
```

**Examples**

```julia
R, (x , y) = PolynomialRing(QQ, ["x", "y"])

I = Ideal(R, x^2 + 1, x*y)

V = lead(I)
```

### Saturation

```@docs
saturation{T <: Nemo.RingElem}(::sideal{T}, ::sideal{T})
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
std(::sideal; ::Bool)
```

```@docs
satstd{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})
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
```

### Reduction

```@docs
reduce(::sideal, ::sideal)
```

```@docs
reduce(::spoly, ::sideal)
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
eliminate(::sideal, ::spoly...)
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
fres{T <: Nemo.RingElem}(::sideal{T}, ::Int, ::String)
```

```@docs
sres{T <: Nemo.RingElem}(::sideal{T}, ::Int)
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
jet(::sideal, ::Int)
```

```@docs
jacobi(::sideal)
```

**Examples**

```julia
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

I = Ideal(R, x^5 - y^2, y^3 - x^6 + z^3)

J1 = jet(I, 3)
J2 = jacobi(I)
```

### Operations on zero-dimensional ideals

```@docs
vdim(::sideal)
```

```@docs
kbase(::sideal)
```

```@docs
highcorner(::sideal)
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

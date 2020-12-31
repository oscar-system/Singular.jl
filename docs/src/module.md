```@meta
CurrentModule = Singular
```

```@meta
DocTestSetup = quote
  using Singular
  using Singular: PolynomialRing, vector, std, reduce
end
```

# Finitely generated modules

Singular.jl allows the creation of submodules of a free module over a Singular polynomial
ring, given by a finite generating set. These are internally stored as a list of elements
of a free module over a polynomial ring $R$. This list of generators can also 
have the property of being a Groebner basis.

The default finitely generated module type in Singular.jl is the Singular `smodule` type.

Module objects have a parent object which represents the class of $R$-modules they belong
to, the data for which is given by the polynomial ring $R$ over which the modules are
defined.

The types of modules and associated parent objects are given in the following table
according to the library providing them.

 Library        | Element type    | Parent type
----------------|-----------------|--------------------------
Singular        | `smodule{T}`    | `Singular.ModuleClass{T}`

These types are parameterised by the type of elements in the polynomial ring $R$.

All module types belong directly to the abstract type `Module{T}` and
all the module class parent object types belong to the abstract type `Set`.

## Module functionality

Singular.jl modules implement standard operations one would expect on modules.
These include:

 * Operations common to all AbstractAlgebra objects, such as `parent`, `base_ring`,
   `elem_type`, `parent_type`, `parent`, `deepcopy`, etc.

Below, we describe all of the functionality for Singular.jl modules that is not included
in this list of basic operations.

### Constructors

Given a Singular polynomial ring $R$, the following constructors are available for
creating modules.

```julia
Module{T <: Nemo.RingElem}(R::PolyRing{T}, vecs::svector{spoly{T}}...)
```

Construct the module over the polynomial ring $R$ whose generators are given 
by the given parameter list of vectors (of length $n$), each component of which is a
polynomial. These vectors represent elements of the free module $R^n$.

Note that `Module` must be prepended with the package name `Singular` to disambiguate
from `Base.Module`.

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

v1 = vector(R, x + 1, x*y + 1, y)
v2 = vector(R, x^2 + 1, 2x + 3y, x)

M = Singular.Module(R, v1, v2)
```

### Basic manipulation

```@docs
ngens(::smodule)
```

```@docs
rank(::smodule)
```

Singular.jl overloads the `setindex!` and `getindex` functions so that one can access
the generators of a module using array notation. Each entry is a vector in $R^n$.

```julia
M[n::Int]
```

```@docs
iszero(::smodule)
```

**Examples**

```
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

v1 = vector(R, x + 1, x*y + 1, y)
v2 = vector(R, x^2 + 1, 2x + 3y, x)

M = Singular.Module(R, v1, v2)

iszero(M) == false
M[1] == v1
n = rank(M)
d = ngens(M)
```

### Standard basis

```@docs
std(::smodule; ::Bool)
```

```@docs
lift_std(::smodule; ::Bool)
```

```@docs
lift_std_syz(::smodule; ::Bool)
```
**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

v1 = vector(R, x + 1, x*y + 1, y)
v2 = vector(R, x^2 + 1, 2x + 3y, x)
v3 = x*v1 + y*v2 + vector(R, x, y + 1, y^2)

M = Singular.Module(R, v1, v2, v3)

G = std(M; complete_reduction=true)
```

### Reduction

```@docs
reduce(::smodule, ::smodule)
```

**Examples**

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Singular Polynomial Ring (QQ),(x,y,z),(dp(3),C), spoly{n_Q}[x, y, z])

julia> v1 = vector(R, R(0), z, -y)
-y*gen(3)+z*gen(2)

julia> v2 = vector(R, -z, R(0), x)
x*gen(3)-z*gen(1)

julia> v3 = vector(R, y, x, R(0))
x*gen(2)+y*gen(1)

julia> v = y*v1+x*v2+z*v3
x^2*gen(3)-y^2*gen(3)+x*z*gen(2)-x*z*gen(1)+y*z*gen(2)+y*z*gen(1)

julia> M = Singular.Module(R, v1, v2, v3)
Singular Module over Singular Polynomial Ring (QQ),(x,y,z),(dp(3),C), with Generators:
-y*gen(3)+z*gen(2)
x*gen(3)-z*gen(1)
x*gen(2)+y*gen(1)

julia> B = std(M; complete_reduction=true)
Singular Module over Singular Polynomial Ring (QQ),(x,y,z),(dp(3),C), with Generators:
y*gen(3)-z*gen(2)
x*gen(2)+y*gen(1)
x*gen(3)-z*gen(1)
y*z*gen(1)

julia> V = Singular.Module(R, v)
Singular Module over Singular Polynomial Ring (QQ),(x,y,z),(dp(3),C), with Generators:
x^2*gen(3)-y^2*gen(3)+x*z*gen(2)-x*z*gen(1)+y*z*gen(2)+y*z*gen(1)

julia> reduce(V,B)
Singular Module over Singular Polynomial Ring (QQ),(x,y,z),(dp(3),C), with Generators:
0

```

### Syzygies

```@docs
syz(::smodule)
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

v1 = vector(R, (x + 1)*y, (x*y + 1)*y, y)
v2 = vector(R, (x + 1)*x, (x*y + 1)*x, x)

M = Singular.Module(R, v1, v2)

Z = syz(M)
```

### Free resolutions

```@docs
sres{T <: Nemo.RingElem}(::smodule{T}, ::Int)
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

v1 = vector(R, x + 1, x*y + 1, y)
v2 = vector(R, x^2 + 1, 2x + 3y, x)

M = std(Singular.Module(R, v1, v2))

F = sres(M, 0)

M1 = Singular.Matrix(M)
M2 = Singular.Matrix(F[2])

# test we have a complex
iszero(M1*M2)
```

### Jet of module

```@docs
jet(::smodule, ::Int)
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"])

v1 = vector(R, x + 1, x*y + 1, y)
v2 = vector(R, x^5 + 1, 2x^3 + 3y^2, x^2)

M = Singular.Module(R, v1, v2)
N = jet(M,3)
```

### Operations over local rings

If the base ring `R` is a local ring, a minimal generating set can be computed
using the following function

```@docs
minimal_generating_set(::smodule)
```

**Examples**

```julia
R, (x, y) = PolynomialRing(QQ, ["x", "y"]; ordering=:negdegrevlex)

has_local_ordering(R) == true

v1 = vector(R, x, y^2)
v2 = vector(R, y - x, y - y^2)
v3 = v1 + v2

M = Singular.Module(R, v1, v2, v3)
   
min = minimal_generating_set(M)
```

```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

```@meta
DocTestSetup = quote
  using Singular
  using Singular: polynomial_ring, vector, std, reduce
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

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular polynomial ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> v1 = vector(R, x + 1, x*y + 1, y)
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)

julia> v2 = vector(R, x^2 + 1, 2x + 3y, x)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)

julia> M = Singular.Module(R, v1, v2)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)
```

### Basic manipulation

```@docs
number_of_generators(::smodule)
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

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular polynomial ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> v1 = vector(R, x + 1, x*y + 1, y)
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)

julia> v2 = vector(R, x^2 + 1, 2x + 3y, x)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)

julia> M = Singular.Module(R, v1, v2)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)

julia> iszero(M) == false
true

julia> M[1] == v1
true

julia> n = rank(M)
3

julia> d = number_of_generators(M)
2
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

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular polynomial ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> v1 = vector(R, x + 1, x*y + 1, y)
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)

julia> v2 = vector(R, x^2 + 1, 2x + 3y, x)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)

julia> v3 = x*v1 + y*v2 + vector(R, x, y + 1, y^2)
x^2*y*gen(2)+x^2*y*gen(1)+x^2*gen(1)+2*x*y*gen(3)+2*x*y*gen(2)+y^2*gen(3)+3*y^2*gen(2)+x*gen(2)+2*x*gen(1)+y*gen(2)+y*gen(1)+gen(2)

julia> M = Singular.Module(R, v1, v2, v3)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)
x^2*y*gen(2)+x^2*y*gen(1)+x^2*gen(1)+2*x*y*gen(3)+2*x*y*gen(2)+y^2*gen(3)+3*y^2*gen(2)+x*gen(2)+2*x*gen(1)+y*gen(2)+y*gen(1)+gen(2)

julia> G = std(M; complete_reduction=true)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
y^2*gen(3)+x*gen(1)+y*gen(2)+gen(2)
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)
```

### Reduction

```@docs
reduce(::smodule, ::smodule)
```

**Examples**

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Singular polynomial ring (QQ),(x,y,z),(dp(3),C), spoly{n_Q}[x, y, z])

julia> v1 = vector(R, R(0), z, -y)
-y*gen(3)+z*gen(2)

julia> v2 = vector(R, -z, R(0), x)
x*gen(3)-z*gen(1)

julia> v3 = vector(R, y, x, R(0))
x*gen(2)+y*gen(1)

julia> v = y*v1+x*v2+z*v3
x^2*gen(3)-y^2*gen(3)+x*z*gen(2)-x*z*gen(1)+y*z*gen(2)+y*z*gen(1)

julia> M = Singular.Module(R, v1, v2, v3)
Singular module over Singular polynomial ring (QQ),(x,y,z),(dp(3),C), with generators:
-y*gen(3)+z*gen(2)
x*gen(3)-z*gen(1)
x*gen(2)+y*gen(1)

julia> B = std(M; complete_reduction=true)
Singular module over Singular polynomial ring (QQ),(x,y,z),(dp(3),C), with generators:
y*gen(3)-z*gen(2)
x*gen(2)+y*gen(1)
x*gen(3)-z*gen(1)
y*z*gen(1)

julia> V = Singular.Module(R, v)
Singular module over Singular polynomial ring (QQ),(x,y,z),(dp(3),C), with generators:
x^2*gen(3)-y^2*gen(3)+x*z*gen(2)-x*z*gen(1)+y*z*gen(2)+y*z*gen(1)

julia> reduce(V,B)
Singular module over Singular polynomial ring (QQ),(x,y,z),(dp(3),C), with generators:
0

```

### Syzygies

```@docs
syz(::smodule)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular polynomial ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> v1 = vector(R, (x + 1)*y, (x*y + 1)*y, y)
x*y^2*gen(2)+x*y*gen(1)+y*gen(3)+y*gen(2)+y*gen(1)

julia> v2 = vector(R, (x + 1)*x, (x*y + 1)*x, x)
x^2*y*gen(2)+x^2*gen(1)+x*gen(3)+x*gen(2)+x*gen(1)

julia> M = Singular.Module(R, v1, v2)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*y^2*gen(2)+x*y*gen(1)+y*gen(3)+y*gen(2)+y*gen(1)
x^2*y*gen(2)+x^2*gen(1)+x*gen(3)+x*gen(2)+x*gen(1)

julia> Z = syz(M)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*gen(1)-y*gen(2)
```

### Free resolutions

```@docs
sres{T <: Singular.FieldElem}(::smodule{spoly{T}}, ::Int)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular polynomial ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> v1 = vector(R, x + 1, x*y + 1, y)
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)

julia> v2 = vector(R, x^2 + 1, 2x + 3y, x)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)

julia> M = std(Singular.Module(R, v1, v2))
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)
x^2*gen(1)+x*gen(3)+2*x*gen(2)+3*y*gen(2)+gen(1)

julia> F = sres(M, 0)
Singular resolution: R^3 <- R^2

julia> M1 = Singular.Matrix(M)
[x + 1, x^2 + 1
x*y + 1, 2*x + 3*y
y, x]

julia> M2 = Singular.Matrix(F[2])
[0
0]

julia> iszero(M1*M2) # test we have a complex
true
```

### Jet of module

```@docs
jet(::smodule, ::Int)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular polynomial ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> v1 = vector(R, x + 1, x*y + 1, y)
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)

julia> v2 = vector(R, x^5 + 1, 2x^3 + 3y^2, x^2)
x^5*gen(1)+2*x^3*gen(2)+x^2*gen(3)+3*y^2*gen(2)+gen(1)

julia> M = Singular.Module(R, v1, v2)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)
x^5*gen(1)+2*x^3*gen(2)+x^2*gen(3)+3*y^2*gen(2)+gen(1)

julia> N = jet(M,3)
Singular module over Singular polynomial ring (QQ),(x,y),(dp(2),C), with generators:
x*y*gen(2)+x*gen(1)+y*gen(3)+gen(2)+gen(1)
2*x^3*gen(2)+x^2*gen(3)+3*y^2*gen(2)+gen(1)
```

### Operations over local rings

If the base ring `R` is a local ring, a minimal generating set can be computed
using the following function

```@docs
minimal_generating_set(::smodule)
```

**Examples**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]; ordering=:negdegrevlex)
(Singular polynomial ring (QQ),(x,y),(ds(2),C), spoly{n_Q}[x, y])

julia> has_local_ordering(R) == true
true

julia> v1 = vector(R, x, y^2)
x*gen(1)+y^2*gen(2)

julia> v2 = vector(R, y - x, y - y^2)
-x*gen(1)+y*gen(2)+y*gen(1)-y^2*gen(2)

julia> v3 = v1 + v2
y*gen(2)+y*gen(1)

julia> M = Singular.Module(R, v1, v2, v3)
Singular module over Singular polynomial ring (QQ),(x,y),(ds(2),C), with generators:
x*gen(1)+y^2*gen(2)
-x*gen(1)+y*gen(2)+y*gen(1)-y^2*gen(2)
y*gen(2)+y*gen(1)

julia> min = minimal_generating_set(M)
2-element Vector{svector{spoly{n_Q}}}:
 y*gen(2)+y*gen(1)
 x*gen(1)-y*gen(2)-y*gen(1)+y^2*gen(2)
```

```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Resolutions

Functions for creating free resolutions of modules and ideals in Singular.jl return a
special Singular object of type `sresolution{T}`. The support in Singular.jl for this
type primarily exists to allow interaction with such resolutions. Free resolutions
can have the property of being minimal, which is specified by the `minimal` field of the
`sresolution{T}` type.

Resolution objects have a parent object which represents the set of resolutions they belong
to, the data for which is given by the polynomial ring $R$ over which the modules in the
resolution are defined.

The types of resolutions and associated parent objects are given in the following table
according to the library providing them.

 Library        | Element type     | Parent type
----------------|------------------|--------------------------
Singular        | `sresolution{T}` | `Singular.ResolutionSet{T}`

These types are parameterised by the type of elements in the polynomial ring $R$ over
which the modules belonging to the resolution are defined.

All resolution types belong directly to the abstract type `SetElem` and
all the resolution set parent object types belong to the abstract type `Set`.

## Resolution functionality

Singular.jl resolutions implement standard operations one would expect on all
AbstractAlgebra compatible objects.
These include:

 * Operations common to all AbstractAlgebra objects, such as `parent`, `base_ring`,
   `elem_type`, `parent_type`, `parent`, `deepcopy`, etc.

Below, we describe all of the functionality for Singular.jl resolutions that is not
included in this list of basic operations.

### Constructors

There are currently two ways to create resolutions in Singular.jl:
They can either be created by taking the free resolution of an ideal or module
over a polynomial ring, as described in the relevant sections of the
documentation, or they can be created by the following constructor:

```@docs
Resolution(::Vector{smodule{T}}) where T <: AbstractAlgebra.FieldElem
```

**Example**

```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Singular Polynomial Ring (QQ),(x,y),(dp(2),C), spoly{n_Q}[x, y])

julia> M1 = Singular.Module(R, vector(R, x), vector(R, y))
Singular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:
x*gen(1)
y*gen(1)

julia> M2 = Singular.Module(R, vector(R, y, -x))
Singular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:
-x*gen(2)+y*gen(1)

julia> F = Resolution([M1, M2]);

julia> F[1]
Singular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:
x*gen(1)
y*gen(1)

julia> F[2]
Singular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:
-x*gen(2)+y*gen(1)
```

Alternatively, resolutions can be refined to minimal resolutions, as described below.

Other than this, there are currently no additional ways to create resolutions in
Singular.jl.

### Basic manipulation

```@docs
length(::sresolution)
```

Singular.jl overloads the `getindex` function so that one can access the modules in a
resolution $F$.

```julia
F[n::Int]
```

**Examples**

```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Singular Polynomial Ring (QQ),(w,x,y,z),(dp(4),C), spoly{n_Q}[w, x, y, z])

julia> I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
Singular ideal over Singular Polynomial Ring (QQ),(w,x,y,z),(dp(4),C) with generators (w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)

julia> F = fres(std(I), 0)
Singular Resolution:
R^1 <- R^5 <- R^6 <- R^2

julia> n = length(F)
3

julia> M1 = F[1]
Singular ideal over Singular Polynomial Ring (QQ),(w,x,y,z),(dp(4),C) with generators (y^2 - w*z, x*y - z^2, x^2 - w*y, w*x - y*z, w^2 - x*z)
```

### Betti numbers

```@docs
betti(::sresolution)
```

**Examples**

```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Singular Polynomial Ring (QQ),(w,x,y,z),(dp(4),C), spoly{n_Q}[w, x, y, z])

julia> I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
Singular ideal over Singular Polynomial Ring (QQ),(w,x,y,z),(dp(4),C) with generators (w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)

julia> F = fres(std(I), 3)
Singular Resolution:
R^1 <- R^5 <- R^6 <- R^2

julia> M = minres(F)
Singular Resolution:
R^1 <- R^5 <- R^5 <- R^1

julia> B = betti(M)
3×4 Matrix{Int32}:
 1  0  0  0
 0  5  5  0
 0  0  0  1
```

### Minimal resolutions

```@docs
minres{T <: Nemo.FieldElem}(::sresolution{spoly{T}})
```

**Examples**

```jldoctest
julia> R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])
(Singular Polynomial Ring (QQ),(w,x,y,z),(dp(4),C), spoly{n_Q}[w, x, y, z])

julia> I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
Singular ideal over Singular Polynomial Ring (QQ),(w,x,y,z),(dp(4),C) with generators (w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)

julia> F = fres(std(I), 3)
Singular Resolution:
R^1 <- R^5 <- R^6 <- R^2

julia> M = minres(F)
Singular Resolution:
R^1 <- R^5 <- R^5 <- R^1
```


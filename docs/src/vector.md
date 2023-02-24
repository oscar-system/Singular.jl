```@meta
CurrentModule = Singular
```

# Free modules and vectors

As generators of finitely generated modules in Singular.jl are given as submodule of
free modules over a polynomial ring $R$, Singular.jl supports creation of the free
module $R^n$ and vectors of length $n$ in such a module.

The Singular.jl type for a vector is `svector{T}`. For the most part, these exist to
help interact with the `smodule{T}` type provided by Singular.

The types of vectors and associated parent objects are given in the following table
according to the library providing them.

 Library        | Element type    | Parent type
----------------|-----------------|--------------------------
Singular        | `svector{T}`    | `Singular.FreeMod{T}`

These types are parameterised by the type of elements in the polynomial ring $R$.

All free module types belong directly to the abstract type `Module{T}` and vector types
belong directly to `ModuleElem{T}`.

## Free module and vector functionality

Singular.jl modules implement standard operations one would expect on vectors and their
associated parent modules.

These include:

 * Operations common to all AbstractAlgebra objects, such as `parent`, `base_ring`,
   `elem_type`, `parent_type`, `parent`, `deepcopy`, etc.

 * Arithmetic operations on vectors: (unary) `-`, `+`, `-`

 * Scalar multiplication of vectors by polynomials, constants and integers

 * Comparison operators: `==`

Below, we describe all of the functionality for Singular.jl free modules that is not
included in this list of basic operations.

### Constructors

Given a Singular polynomial ring $R$ and a rank $n$, the following constructors are
available for creating free modules.

```julia
FreeModule{T <: AbstractAlgebra.RingElem}(R::PolyRing{T}, n::Int)
```

Construct the free module $R^n$ over the polynomial ring $R$. Elements of the module
returned by this function are vectors of length $n$, with entries in $R$.

```julia
vector{T <: AbstractAlgebra.RingElem}(R::PolyRing{T}, coords::spoly{T}...)
```

Construct the vector whose coordinates (which are elements of $R$) are the given
parameters.

**Examples**

```julia
R, (x, y) = polynomial_ring(QQ, ["x", "y"])

M = FreeModule(R, 3)
v2 = M([x + 1, x*y + 1, y])

v1 = vector(R, x + 1, x*y + 1, y)
```

### Basic manipulation


```@docs
rank(::FreeMod)
```

```@docs
gens{T <: AbstractAlgebra.RingElem}(::FreeMod{T})
```

**Examples**

```
R, (x, y) = polynomial_ring(QQ, ["x", "y"])

M = FreeModule(R, 5)

v = gens(M)
r = rank(M)
```

### Conversions

To convert the internal Singular representation of an `svector{T}` to a Julia array
whose entries have type `T`, we have the following conversion routine.

```julia
Array{T <: Nemo.RingElem}(v::svector{spoly{T}})
```

**Examples**

```julia
R, (x, y) = polynomial_ring(QQ, ["x", "y"])

v1 = vector(R, x + 1, x*y + 1, y)

V = Array(v1)
```

### Jet of vectors

```@docs
jet{T <: AbstractAlgebra.RingElem}(::svector{spoly{T}}, ::Int)
```

**Examples**
```julia
R, (x, y) = polynomial_ring(QQ, ["x", "y"])

v = vector(R, x^5 + 1, 2x^3 + 3y^2, x^2)
w = jet(v, 3)
```


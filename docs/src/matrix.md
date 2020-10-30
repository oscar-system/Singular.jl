```@meta
CurrentModule = Singular
```

# Matrices

Singular internally allows for matrices over polynomial rings to be created extremely
efficiently from ideals and modules (often without copying data). This allows for
introspection of modules and operations that can be expressed in terms of matrices (e.g.
composition of $R$-module homomorphisms) to be computed, at a low level.

The default matrix type in Singular.jl is the `smatrix` type.

Matrix objects have a parent object which represents the space of matrices they belong
to, the data for which is given by the polynomial ring $R$ over which the matrices are
defined, and the number of rows and columns of the matrices in the space.

The types of matrices and associated parent objects are given in the following table
according to the library providing them.

 Library        | Element type    | Parent type
----------------|-----------------|--------------------------
Singular        | `smatrix{T}`    | `Singular.MatrixSpace{T}`

These types are parameterised by the type of elements in the polynomial ring $R$ over
which the matrices are defined.

All matrix types belong directly to the abstract type `SetElem` and
all the matrix space parent object types belong to the abstract type `Set`.

## Matrix functionality

Singular.jl matrices implement standard operations one would expect.
These include:

 * Operations common to all AbstractAlgebra objects, such as `parent`, `base_ring`,
   `elem_type`, `parent_type`, `parent`, `deepcopy`, etc.

The following parts of the Matrix interface from AbstractAlgebra are also implemented:

  * arithmetic operations: `+`, `-`, `*`
  * comparison: `==`

## Basic constructors

Using Singular.jl, the user has two basic constructors for matrices.

```@docs
identity_matrix(::PolyRing, ::Int)
```

```@docs
zero_matrix(::PolyRing, ::Int, ::Int)
```

**Examples**
```julia
R, (x, y, u, v, w) = Singular.PolynomialRing(Singular.QQ, ["x", "y", "u", "v", "w"])

identity_matrix(R, 4)

zero_matrix(R, 3, 8)
```

## Basic manipulation

The following functions are provided to manipulate matrices and to set and
retrieve entries other basic data associated with the matrices.

```@docs
nrows(::smatrix)
```

```@docs
ncols(::smatrix)
```

```@docs
getindex(::smatrix{T}, ::Int, ::Int) where T <: AbstractAlgebra.RingElem
```

```@docs
setindex!(::smatrix, ::spoly, ::Int, ::Int)
```

```@docs
transpose(M::smatrix{T}) where T <: AbstractAlgebra.RingElem
```

```@docs
iszero(M::smatrix)
```

**Examples**
```julia
R, (x, y, u, v, w) = Singular.PolynomialRing(QQ, ["x", "y", "u", "v", "w"])

M = identity_matrix(R, 4)

nrows(M)

ncols(M)

iszero(M)

M[3, 4] = x*y + 5*u*w

N = transpose(M)

N[4, 3]
```

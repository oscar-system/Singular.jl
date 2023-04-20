```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Nemo rings and fields

Any type that satisfies AbstractAlgebra.jl Ring or Field interface, such as all Nemo
ring and field types, can be used as coefficient rings in Singular.jl. These are
implemented via the Singular `n_RingElem{T}` and `n_FieldElem{T}` types, parameterised
by the given Nemo/AbstractAlgebra element type.

The associated parent object of type `N_Ring{T}` or `N_Field{T}` can be constructed by
a call to the `CoefficientRing` constructor. In practice, however, this constructor is
only used internally, and Nemo rings and fields work directly as Singular coefficient
rings, and all the coercions and ad hoc functions that one would expect to be present
are implemented.

The Singular.jl `n_RingElem` (resp. `n_FieldElem`) types belong directly to the
abstract type `RingElem` (resp. `FieldElem`) and their parent object types belong
to the abstract type `Ring` (resp. `Field`). We also have the following type
definitions for compatibility with previous version of Singular.jl.

```julia
const N_unknown{T} = Union{N_Ring{T}, N_Field{T}}
const n_unknown{T} = Union{n_RingElem{T}, n_FieldElem{T}}
```

All of the Singular polynomial arithmetic should work for any Nemo ring and everything,
including ideals, modules, standard basis, syzygies, resolutions, etc., should work
with any Nemo field.

Specialised efficient wrappers exist for certain Nemo coefficient ring types.

## Nemo ring functionality

Singular.jl foreign ring types provide all of the AbstractAlgebra defined ring and
some field functionality.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/ring>

<https://nemocas.github.io/AbstractAlgebra.jl/latest/field>

Parts of the Euclidean Ring interface may also be implemented, though Singular will
report an error if division is meaningless (even after cancelling zero divisors).

<https://nemocas.github.io/AbstractAlgebra.jl/latest/euclidean_interface>

Below, we describe the functionality that is specific to the Singular foreign ring
interface that is not already listed at the given links.

### Constructors

Given an AbstractAlgebra compatible ring $R$, e.g. a Nemo ring, we have the following
constructor, which returns the associated Singular.jl coefficient ring.

```julia
CoefficientRing(R::Ring)
```

If there are generators to be coerced from Nemo/AbstractAlgebra into corresponding
elements, the Singular.jl coefficient ring can be used to coerce them to a Singular
`n_RingElem` or `n_FieldElem` element.

**Examples**

```julia
R, x = Nemo.polynomial_ring(ZZ, "x")
S = CoefficientRing(R)
t = S(x)
```

Note that it is unlikely that a user directly needs to construct the Singular
coefficient ring from a Nemo ring, since the Singular.jl constructors are designed to
accept Nemo coefficient rings directly. Singular.jl automatically constructs the
required Singular coefficient ring and makes use of it.


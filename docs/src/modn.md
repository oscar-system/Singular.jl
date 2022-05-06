```@meta
CurrentModule = Singular
```

# Integers mod n

Integers mod $n$ are implemented via the Singular `n_Zn` type for any positive modulus
that can fit in a Julia `Int`.

The associated ring of integers mod $n$ is represented by a parent object which can
be constructed by a call to the `ResidueRing` constructor.

The types of the parent objects and elements of the associated rings of integers modulo
n are given in the following table according to the library providing them.

 Library        | Element type  | Parent type
----------------|---------------|--------------------
Singular        | `n_Zn`        | `Singular.N_ZnRing`

All integer mod $n$ element types belong directly to the abstract type `RingElem` and
all the parent object types belong to the abstract type `Ring`.

## Integer mod $n$ functionality

Singular.jl integers modulo $n$ provide all the AbstractAlgebra ring and residue ring
functionality.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/ring>

<https://nemocas.github.io/AbstractAlgebra.jl/latest/residue>

Parts of the Euclidean Ring interface may also be implemented, though Singular will
report an error if division is meaningless (even after cancelling zero divisors).

<https://nemocas.github.io/AbstractAlgebra.jl/latest/euclidean_interface>

Below, we describe the functionality that is specific to the Singular integers mod $n$
ring and not already listed at the given links.

### Constructors

Given a ring $R$ of integers modulo $n$, we also have the following coercions in
addition to the standard ones expected.

```julia
R(n::n_Z)
R(n::fmpz)
```

Coerce a Singular or Flint integer value into the ring.

### Basic manipulation

**Examples**

```
R = ResidueRing(ZZ, 26)
a = R(5)

is_unit(a)
c = characteristic(R)
```


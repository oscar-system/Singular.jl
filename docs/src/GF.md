```@meta
CurrentModule = Singular
```

# Finite fields

Finite fields are implemented via the Singular `n_GF` type for any characteristic and
degree contained in the Singular Conway tables.

The associated finite field is represented by a parent object which can be constructed
by a call to the `FiniteField` constructor.

The types of the parent objects and elements of the associated finite fields are given
in the following table according to the library providing them.

 Library        | Element type  | Parent type
----------------|---------------|--------------------
Singular        | `n_GF`        | `Singular.N_GField`

All finite field element types belong directly to the abstract type `FieldElem` and
all the parent object types belong to the abstract type `Field`.

## Finite field functionality

Singular.jl finite fields implement the Field interface of AbstractAlgebra.jl.

<https://nemocas.github.io/AbstractAlgebra.jl/fields.html>

Below, we describe the functionality that is specific to Singular finite field and not
already listed at the given link.

### Constructors

The following constructors are available to create finite fields and their elements.

```@docs
Singular.FiniteField(::Int, ::Int, ::String; ::Bool)
```

Given a finite field $R$, we also have the following coercions in addition to the
standard ones expected.

```julia
R(n::fmpz)
```

Coerce a Flint integer value into the field.

### Basic manipulation

```@docs
Singular.degree(::N_GField)
```

**Examples**

```julia
R,w = FiniteField(7, 2, "w")
w^48 == 1
a = R(5)

is_unit(a)
c = characteristic(R)
d = degree(R)
```


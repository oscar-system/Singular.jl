```@meta
CurrentModule = Singular
```

# Integers mod p

Integers mod a prime $p$ are implemented via the Singular `n_Zp` type for any positive
prime modulus less than $2^{29}$.

The associated field of integers mod $p$ is represented by a parent object which can
be constructed by a call to the `Fp` constructor.

The types of the parent objects and elements of the associated fields of integers modulo
$p$ are given in the following table according to the library providing them.

 Library        | Element type  | Parent type
----------------|---------------|--------------------
Singular        | `n_Zp`        | `Singular.N_ZpField`

All integer mod $p$ element types belong directly to the abstract type `FieldElem` and
all the parent object types belong to the abstract type `Field`.

## Integer mod $p$ functionality

Singular.jl integers modulo $p$ provides the field and residue ring functionality of
AbstractAlgebra.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/field>

<https://nemocas.github.io/AbstractAlgebra.jl/latest/residue>

Below, we describe the functionality that is specific to the Singular integers mod $p$
field and not already listed at the given links.

### Constructors

The following constructors are available to create the field of integers modulo a
prime $p$.

```julia
Fp(p::Int; cached=true)
```

Construct the field of integers modulo $p$. By default, the field is cached, so that
all fields of integers modulo $p$ have the same parent object. If this is not the
desired behaviour, the `cached` parameter can be set to `false`. If $p$ is not a prime
or $p$ is not in the range $(0, 2^{29})$, an exception is raised.

Given a field $R$ of integers modulo $p$, we also have the following coercions in
addition to the standard ones expected.

```julia
R(n::n_Z)
R(n::fmpz)
```

Coerce a Singular or Flint integer value into the field.

### Basic manipulation

```@docs
isunit(::n_Zp)
```

```@docs
Singular.characteristic(::N_ZpField)
```

**Examples**

```julia
R = Fp(23)
a = R(5)

isunit(a)
c = characteristic(R)
```

### Conversions

```
Int(n::n_Zp)
```

Lift the integer $n$ modulo $p$ to a Julia `Int`. The result is always in the range
$[0, p)$.

**Examples**

```
R = Fp(23)
a = R(5)

b = Int(a)
```


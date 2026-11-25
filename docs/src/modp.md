```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Integers mod p

Integers mod a prime $p$ are implemented via Singular's machine-integer `n_Zp`
representation when the modulus fits and via its GMP-backed modular representation
otherwise.

The associated field of integers mod $p$ is represented by a parent object which can
be constructed by a call to the `Fp` constructor.

The types of the parent objects and elements of the associated fields of integers modulo
$p$ are given in the following table according to the library providing them.

 Library        | Element type  | Parent type
----------------|---------------|--------------------
Singular        | `n_Zp`        | `Singular.N_ZpField`
Singular (GMP)  | `n_Zn`        | `Singular.N_ZnRing`

The GMP-backed representation is shared with general residue rings. Singular marks
its coefficient domain as a field when the modulus is prime.

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
Fp(p::Integer; cached=true)
```

Construct the field of integers modulo $p$. Singular chooses the coefficient
representation according to the size of $p$. By default, the result is cached, so
that repeated constructions use the same parent object. If this is not desired, the
`cached` parameter can be set to `false`. If $p$ is not positive and prime, an
exception is raised.

Given a field $R$ of integers modulo $p$, we also have the following coercions in
addition to the standard ones expected.

```julia
R(n::n_Z)
R(n::ZZRingElem)
```

Coerce a Singular or Flint integer value into the field.

### Basic manipulation

**Examples**

```jldoctest
julia> R = Fp(23)
Finite field of characteristic 23

julia> a = R(5)
5

julia> is_unit(a)
true

julia> c = characteristic(R)
23
```

### Conversions

```
Int(n::n_Zp)
```

Lift the integer $n$ modulo $p$ to a Julia `Int`. The result is always in the range
$[0, p)$.

**Examples**

```jldoctest
julia> R = Fp(23)
Finite field of characteristic 23

julia> a = R(5)
5

julia> b = Int(a)
5
```

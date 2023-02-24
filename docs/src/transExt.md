```@meta
CurrentModule = Singular
```

# Function fields

Function fields are implemented via the Singular `n_transExt` type for prime fields of any characteristic.

The associated function field is represented by a parent object which can be constructed
by a call to the `FunctionField` constructor.

The types of the parent objects and elements of the associated function fields are given
in the following table according to the library providing them.

 Library        | Element type  | Parent type
----------------|---------------|--------------------
Singular        | `n_transExt`        | `Singular.N_FField`

All function field element types belong directly to the abstract type `FieldElem` and
all the parent object types belong to the abstract type `Field`.

## Function field functionality

Singular.jl function fields provide all the functionality for fields described by
AbstractAlgebra.jl.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/field>

Below, we describe the functionality that is specific to Singular function field and not
already listed at the given link.

### Constructors

The following constructors are available to create function fields and their elements.

```@docs
Singular.FunctionField(::Field, ::Vector{String}; cached::Bool = true)
Singular.FunctionField(::Field, ::Vector{Symbol}; cached::Bool = true)
```

In case the user does not want to specify a transcendence basis the following
constructor can be used.

```@docs
Singular.FunctionField(::Field, ::Int; cached::Bool = true)
```

Given a function field $F$, we also have the following coercions in addition to the
standard ones expected.

```julia
F(n::ZZRingElem)
```

Coerce a Flint integer value into the field.

**Examples**

```julia
F1, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])

x1 = a*b + c

F2, (a1, a2, a3) = FunctionField(Fp(5), 3)

x2 = a1^5 + a2*a3^4
```

### Basic manipulation

```@docs
numerator(::n_transExt)
```

```@docs
denominator(::n_transExt)
```

```@docs
transcendence_degree(::N_FField)
```

```@docs
transcendence_basis(::N_FField)
```

```@docs
n_transExt_to_spoly(x::n_transExt; parent::PolyRing)
```

**Examples**

```julia
F1, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])
x = F1(5)*a
y = a^2 *b+a*b+b^2

is_unit(x)
char = characteristic(F1)
d = transcendence_degree(F1)

S, = polynomial_ring(QQ, ["a", "b", "c"])

p = n_transExt_to_spoly(y, parent_ring = S)

F2, = FunctionField(Fp(7), 4)
B = transcendence_basis(F2)
```


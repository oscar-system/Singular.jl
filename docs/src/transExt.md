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

Singular.jl function fields implement the Field interface of AbstractAlgebra.jl.

[https://nemocas.github.io/AbstractAlgebra.jl/fields.html](https://nemocas.github.io/AbstractAlgebra.jl/fields.html)

Below, we describe the functionality that is specific to Singular function field and not
already listed at the given link.

### Constructors

The following constructors are available to create function fields and their elements.

```@docs
Singular.FunctionField(::Field, ::Array{String, 1}; ::Bool)
```

In case the user does not want to specify a trancendece basis the following 
constructor can be used.

```@docs
Singular.FunctionField(::Field, ::Int; ::Bool)
```

Given a function field $F$, we also have the following coercions in addition to the
standard ones expected.

```julia
F(n::fmpz)
```

Coerce a Flint integer value into the field.

### Basic manipulation

```@docs
numerator(::n_transExt)
```

```@docs
denominator(::n_transExt)
```

```@docs
Singular.transcendence_degree(::N_FField)
```

```@docs
Singular.transcendence_basis(::N_FField)
```

```@docs
Singular.characteristic(::N_FField)
```

```@docs
isunit(::n_transExt)
```

**Examples**

```julia
F1, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])
x = F(5)

isunit(x)
char = characteristic(F1)
d = transcendence_degree(F1)

F2, = FunctionField(Fp(7), 4)
B = transcendence_basis(F2)
```


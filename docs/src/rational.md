```@meta
CurrentModule = Singular
```

# Rational field

Singular.jl provides rational numbers via Singular's `n_Q` type.

There is a constant parent object representing the field of rationals, called `QQ`
in Singular.jl. It is defined by `QQ = Rationals()`, which calls the constructor for
the unique field of rationals in Singular.

## Rational functionality

The rationals in Singular.jl implement the Field interface defined by AbstractAlgebra.jl.
They also implement the Fraction Field interface.

[https://nemocas.github.io/AbstractAlgebra.jl/fields.html](https://nemocas.github.io/AbstractAlgebra.jl/fields.html)

[https://nemocas.github.io/AbstractAlgebra.jl/fraction_fields.html](https://nemocas.github.io/AbstractAlgebra.jl/fraction_fields.html)

We describe here only the extra functionality provided by Singular that is not already
described in those interfaces.

### Constructors

In addition to the standard constructors required for the interfaces listed above,
Singular.jl provides the following constructors.

```
QQ(n::n_Z)
QQ(n::fmpz)
```

Construct a Singular rational from the given integer $n$.

### Basic manipulation

```@docs
numerator(::n_Q)
```

```@docs
denominator(::n_Q)
```

```@docs
abs(::n_Q)
```

```julia
f = QQ(-12, 7)

h = numerator(QQ)
k = denominator(QQ)
m = abs(f)
```

### Comparison

Here is a list of the comparison functions implemented, with the understanding
that `isless` provides all the usual comparison operators.

Function                   |
---------------------------|
`isless(a::n_Q, b::n_Q)`   |

We also provide the following ad hoc comparisons which again provide all of the
comparison operators mentioned above.

Function                     |
-----------------------------|
`isless(a::n_Q, b::Integer)` |
`isless(a::Integer, b::n_Q)` |

**Examples**

```julia
a = QQ(12, 7)
b = QQ(-3, 5)

a > b
a != b
a > 1
5 >= b
```

### Rational reconstruction

```@docs
reconstruct(::n_Z, ::n_Z)
```

The following ad hoc versions of the same function also exist.

```julia
reconstruct(::n_Z, ::Integer)
reconstruct(::Integer, ::n_Z)
```

**Examples**

```julia
q1 = reconstruct(ZZ(7), ZZ(3))
q2 = reconstruct(ZZ(7), 5)
```

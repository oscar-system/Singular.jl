```@meta
CurrentModule = Singular
```

# Integers

The default integer type in Singular.jl is the Singular `n_Z` integer type.

The associated ring of integers is represented by the constant parent object which can
be constructed by a call to `Singular.Integers()`.

For convenience we define

```
ZZ = Singular.Integers()
```

so that integers can be constructed using `ZZ`. Note that this is the name of a
specific parent object, not the name of its type.

The types of the integer ring parent objects and elements of the associated
rings of integers are given in the following table according to the library
providing them.

 Library        | Element type  | Parent type
----------------|---------------|--------------------
Singular        | `n_Z`         | `Singular.Integers`

All integer element types belong directly to the abstract type `RingElem` and
all the integer ring parent object types belong to the abstract type `Ring`.

## Integer functionality

Singular.jl provides all the ring and possibly some parts of the Euclidean ring
functionality of AbstractAlgebra.jl.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/ring>

<https://nemocas.github.io/AbstractAlgebra.jl/latest/euclidean_interface>

Below, we describe the functionality that is specific to the Singular integer ring.

### Constructors

```julia
ZZ(n::Integer)
```

Coerce a Julia integer value into the integer ring.

### Basic manipulation

```@docs
denominator(::n_Z)
```

```@docs
numerator(::n_Z)
```

**Examples**

```
a = ZZ(-12)

is_unit(a)
n = numerator(a)
d = denominator(a)
c = abs(a)
```

### Euclidean division

Singular.jl provides a number of Euclidean division operations. Recall that
for a dividend $a$ and divisor $b$, we can write $a = bq + r$ with
$0 \leq |r| < |b|$. We call $q$ the quotient and $r$ the remainder.

In the following table we list the division functions and their rounding
behaviour. We also give the return value of the function, with $q$ representing
return of the quotient and $r$ representing return of the remainder.

Function                    | Return | Rounding
----------------------------|--------|------------------------
`divrem(a::n_Z, b::n_Z)`    | q, r   | towards zero
`rem(a::n_Z, b::n_Z)`       | r      | towards zero
`mod(a::n_Z, b::n_Z)`       | r      | down

**Examples**

```julia
a = ZZ(-12)
b = ZZ(5)

q, r = divrem(a, b)
r = mod(a, b)
c = a % b
```

### Comparison

Here is a list of the comparison functions implemented, with the understanding
that `isless` provides all the usual comparison operators.

Function                   |
---------------------------|
`isless(a::n_Z, b::n_Z)`   |

We also provide the following ad hoc comparisons which again provide all of the
comparison operators mentioned above.

Function                     |
-----------------------------|
`isless(a::n_Z, b::Integer)` |
`isless(a::Integer, b::n_Z)` |

**Examples**

```julia
a = ZZ(12)
b = ZZ(3)

a < b
a != b
a > 4
5 <= b
```


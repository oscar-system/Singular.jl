```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Rational field

Singular.jl provides rational numbers via Singular's `n_Q` type.

There is a constant parent object representing the field of rationals, called `QQ`
in Singular.jl. It is defined by `QQ = Rationals()`, which calls the constructor for
the unique field of rationals in Singular.

## Rational functionality

The rationals in Singular.jl provide all functionality for fields and fraction fields
described by AbstractAlgebra.jl.

<https://nemocas.github.io/AbstractAlgebra.jl/latest/field>

<https://nemocas.github.io/AbstractAlgebra.jl/latest/fraction>

We describe here only the extra functionality provided by Singular that is not already
described in those interfaces.

### Constructors

In addition to the standard constructors required for the interfaces listed above,
Singular.jl provides the following constructors.

```
QQ(n::n_Z)
QQ(n::ZZRingElem)
```

Construct a Singular rational from the given integer $n$.

**Examples**

```jldoctest
julia> f = QQ(-12, 7)
-12//7

julia> h = numerator(f)
-12

julia> k = denominator(f)
7

julia> m = abs(f)
12//7

julia> a = QQ(12, 7)
12//7

julia> b = QQ(-3, 5)
-3//5

julia> a > b
true

julia> a != b
true

julia> a > 1
true

julia> 5 >= b
true
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

```jldoctest
julia> q1 = reconstruct(ZZ(7), ZZ(3))
1

julia> q2 = reconstruct(ZZ(7), 5)
-1//2
```

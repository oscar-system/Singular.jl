```@meta
CurrentModule = Singular
```

# Multivariate polynomials

Singular.jl allows the creation of multivariate polynomials over any of the coefficient
rings described above.

The default multivariate polynomial type in Singular.jl is the Singular `spoly` type.

The associated polynomial ring is represented by a parent object which can be
constructed by a call to the `PolynomialRing` constructor.

The types of the polynomial ring parent objects and elements thereof are given in the
following table according to the library provding them.

 Library        | Element type  | Parent type
----------------|---------------|-----------------------
Singular        | `spoly{T}`    | `Singular.PolyRing{T}`

These types are parameterised by the type of elements in the coefficient ring of the
polynomials.

All polynomial types belong directly to the abstract type `RingElem` and
all the polynomial ring parent object types belong to the abstract type `Ring`.

## Multivariate polynomial functionality

Singular.jl polynomials implement the Multivariate Polynomial Ring interface of
AbstractAlgebra.jl.

[https://nemocas.github.io/AbstractAlgebra.jl/mpolynomial_rings.html](https://nemocas.github.io/AbstractAlgebra.jl/mpolynomial_rings.html)

In particular, Singular polynomials are sparse distributed, but do not have random
access. Instead, they implement iterator access to terms. This is due to their storage
in a linked list, for efficient implementation of Groebner basis algorithms.

Some polynomial rings may also implement part of the Euclidean Ring interface, where
this is appropriate.

[https://nemocas.github.io/AbstractAlgebra.jl/euclidean.html](https://nemocas.github.io/AbstractAlgebra.jl/euclidean.html)

Below, we describe the functionality that is specific to the Singular multivariate 
polynomials that is not documented in the general multivariate interface.

### Constructors

```julia
PolynomialRing(R::Union{Ring, Field}, s::Array{String, 1};
   cached::Bool = true, ordering::Symbol = :degrevlex,
      ordering2::Symbol = :comp1min, degree_bound::Int = 0)
```

Returns a tuple, $S, x$ consisting of a multivariate polynomial ring $S$ and an array
of variables (from which polynomials can be constructed). The ring $R$ must be a valid
Singular coefficient ring, or any Nemo/AbstractAlgebra coefficient ring. The array $s$
must be a list of strings corresponding to how the variables will be printed. By default,
there will only be one Singular polynomial ring in the system for each combination of
coefficient ring, list of variable names, ordering and degree bound. This is
accomplished by making use of a global cache. If this is not the desired behaviour,
`false` can be passed to the optional argument `cached`. 

Two orderings can be specified, one for term ordering of the polynomials, and another
for ordering of module components. They can occur in either order, the first taking
precedence over the other, when the polynomials are used to represent module generators.
If either is not specified, the indicated default is used.

The options for polynomial term ordering are the symbols, `:lex`, `:deglex` and
`:degrevlex`, and the options for module component ordering are `comp1min` and
`comp1max`.

If one has an a priori bound on the degree in each variable of a polynomial (including
for all intermediate computations in this ring), one can specify it using the
`degree_bound` optional parameter. Singular may then be able to use a more efficient
representation internally, which will use less memory and allow for faster arithmetic.
By default, Singular uses a bound of 16 bits internally for the exponent of each
variable, however this is a signed value, so that the default is for nonnegative
exponents that fit in 15 bits.

Note that internally, Singular may use a higher bound than specified, if it will not
increase the amount of storage required.

**Examples**

```julia
R, (x, y, z) = PolynomialRing(ZZ, ["x", "y", "z"])

S, vars = PolynomialRing(QQ, ["x", "y"]; ordering=:deglex)

T, x = PolynomialRing(ZZ, ["x$i" for i in 1:5];
       ordering=:comp1max, ordering2=:degrevlex, degree_bound=5)
```

See also the convenience macros below for simple use cases.

The following function allows creating a Singular polynomial ring from a given
polynomial ring of type AbstractAlgebra.Generic.MPolyRing:

```@docs
PolynomialRing(R::AbstractAlgebra.Generic.MPolyRing{T}; cached::Bool = true, ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min, degree_bound::Int = 0)  where {T <: RingElement}
```

### Polynomial ring macros

For convenience, we provide some macros for constructing polynomial rings and injecting
the variables into scope. These are easier to use, but have some limitations. In
particular, they can only be used at the top level by the user, and cannot be used
programmatically or in library code (it is not possible to inject an arbitrary number
of variables into scope inside a function).

The macros are designed for simple use cases, and do not offer the full power of the
most general constructor above.

```julia
@PolynomialRing(R, s, n, o)
```

Given a coefficient ring $R$, a root variable name, e.g. `"x"`, a number of variable
$n$ and a polynomial term ordering `o`, create the variables `x1, x2, ..., xn` and
inject them into scope, and return the corresponding polynomial ring `S`.

```julia
@PolynomialRing(R, s, n)
```

As per the previous macro, with a default of `:degrevlex` for the polynomial term
ordering.

**Examples**

```julia
S = @PolynomialRing(ZZ, "x", 5, :deglex)

T = @PolynomialRing(QQ, "y", 10)
```

### Basic manipulation

```@docs
ngens(::PolyRing)
```

```@docs
has_global_ordering(R::PolyRing)
```

```@docs
Singular.characteristic(R::PolyRing)
```

```@docs
degree_bound(R::PolyRing)
```

```@docs
lead_exponent(p::spoly)
```

**Examples**

```
R = @PolynomialRing(ZZ, "x", 3)

n = ngens(R)
has_global_ordering(R) == true
c = characteristic(R)
L = degree_bound(R)
exps = lead_exponent(x1*x2 + 3x1*x2^2 + x3 + 2)
```

### Content and primitive part

When coefficient rings have a meaningful GCD function, the following functions are
available.

```@docs
Singular.primpart(x::spoly)
```

```@docs
Singular.content(x::spoly)
```

**Examples**

```julia
R = @PolynomialRing(ZZ, "x", 2)

f = 3x1^2 + 3x1*x2 + 6x2^2

p = primpart(f)
c = content(f)
```




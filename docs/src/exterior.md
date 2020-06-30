```@meta
CurrentModule = Singular
```

# Exterior algebras

Singular.jl allows the creation of Exterior algebras over any of the coefficient
rings described above.

The default Exterior algebra type in Singular.jl is the Singular `pexterior` type.

The associated Exterior algebra is represented by a parent object which can be
constructed by a call to the `ExteriorAlgebra` constructor.

The types of the Exterior algebra parent objects and elements thereof are given in the
following table according to the library providing them.

 Library        | Element type      | Parent type
----------------|-------------------|--------------------------
Singular        | `pexterior{T}`    | `Singular.ExteriorAlgebra{T}`

These types are parameterised by the type of elements in the coefficient ring of the
Exterior algebra.

All Exterior algebra element types belong directly to the abstract type
`AbstractAlgebra.NCRingElem` and all the Exterior algebra parent object types belong
to the abstract type `AbstractAlgebra.NCRing`.

### Constructors

```julia
ExteriorAlgebra(R::Union{Ring, Field}, s::Array{String, 1};
   cached::Bool = true, ordering::Symbol = :degrevlex,
      ordering2::Symbol = :comp1min, degree_bound::Int = 0)
```

Returns a tuple, $S, x$ consisting of an Exterior algebra $S$ and an array $x$
of variables (from which elements of the algebra can be constructed).

The ring $R$ must be a valid Singular coefficient ring, or any
Nemo/AbstractAlgebra coefficient ring. The array $s$ must be a list of strings
corresponding to how the variables will be printed. By default,
there will only be one Singular Exterior algebra in the system for each combination of
coefficient ring, list of variable names, ordering and degree bound. This is
accomplished by making use of a global cache. If this is not the desired behaviour,
`false` can be passed to the optional argument `cached`. 

Two orderings can be specified, one for term ordering, and another for ordering
of module components. They can occur in either order, the first taking
precedence over the other, when the algebra elements are used to represent module
generators. If either is not specified, the indicated default is used.

The options for term ordering are the symbols, `:lex`, `:deglex`, `:degrevlex`,
`:neglex`, `:negdeglex` and `:negdegrevlex`, and the options for module component
ordering are `comp1min` and `comp1max`.

If one has an a priori bound on the degree in each variable (including
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
R, (x, y) = ExteriorAlgebra(ZZ, ["x", "y"])

S, vars = ExteriorAlgebra(QQ, ["x", "y"]; ordering=:deglex)

T, x = ExteriorAlgebra(ZZ, ["x$i" for i in 1:5];
       ordering=:comp1max, ordering2=:degrevlex, degree_bound=5)
```

See also the convenience macros below for simple use cases.

### Exterior algebra macros

For convenience, we provide some macros for constructing Exterior algebras and injecting
the variables into scope. These are easier to use, but have some limitations. In
particular, they can only be used at the top level by the user, and cannot be used
programmatically or in library code (it is not possible to inject an arbitrary number
of variables into scope inside a function).

The macros are designed for simple use cases, and do not offer the full power of the
most general constructor above.

```julia
@ExteriorAlgebra(R, s, n, o)
```

Given a coefficient ring $R$, a root variable name, e.g. `"x"`, a number of variable
$n$ and a term ordering `o`, create the variables `x1, x2, ..., xn` and inject them
into scope, and return the corresponding Exterior algebra `S`.

```julia
@ExteriorAlgebra(R, s, n)
```

As per the previous macro, with a default of `:degrevlex` for the term
ordering.

**Examples**

```julia
S = @ExteriorAlgebra(ZZ, "x", 5, :deglex)

T = @ExteriorAlgebra(QQ, "y", 10)
```


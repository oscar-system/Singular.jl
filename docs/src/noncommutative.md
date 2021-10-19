```@meta
CurrentModule = Singular
```

# Noncommutative algebras

Singular.jl allows the creation of various noncommutative algebras over any of
the coefficient rings described above. The constructors of the parent objects
and elements thereof are given in the following table.

 Constructor     | Element type    | Parent type        | SINGULAR kernel subsystem
-----------------|-----------------|--------------------|-----------
GAlgebra         | `sgpoly{T}`     | `GPolyRing{T}`     | PLURAL
WeylAlgebra      | `sweylpoly{T}`  | `WeylPolyRing{T}`  | PLURAL
ExteriorAlgebra  | `sextpoly{T}`   | `ExtPolyRing{T}`   | PLURAL + quotient
FreeAlgebra      | `slppoly{T}`    | `LPPolyRing{T}`    | LETTERPLACE

These types are parameterized by the type of elements in the coefficient ring
of the algebra. All noncommutative algebra element types belong directly to the
abstract type `AbstractAlgebra.NCRingElem` and all the noncommuative algebra
parent object types belong to the abstract type `AbstractAlgebra.NCRing`.

## Constructors

All constructors returns a tuple, $R, x$ consisting of a parent object $R$ and
an array $x$ of variables from which elements of the algebra can be constructed.

For constructors taking an ordering, two orderings can be specified, one for
term ordering, and another for ordering of module components. They can occur in
either order, the first taking precedence over the other, when the algebra
elements are used to represent module generators. If either is not specified,
the indicated default is used. The options for term ordering are the symbols
`:lex`, `:deglex`, `:degrevlex`,`:neglex`, `:negdeglex` and `:negdegrevlex`,
and the options for module component ordering are `comp1min` and `comp1max`.

By default there will only be one parent object in the system for each
combination of arguments. This is accomplished by making use of a global cache.
If this is not the desired behaviour `cached = false` may be passed.

### GAlgebra

```julia
GAlgebra(R::PolyRing{T}, C::smatrix{spoly{T}}, D::smatrix{spoly{T}};
         cached::Bool = true) where T <: Nemo.RingElem
```

Construct the G-algebra from a commutative polynomial ring $R$ and matrices $C$,
$D$ over $R$. If the variables of $R$ are $x_1,\dots,x_n$ then the noncommutative
algebra is constructed with relations $x_j x_i = c_{i,j} x_i x_j + d_{i,j}$ for
$1 \le i < j \le n$. The $c_{i,j}$ must be constant polynomials and the monomial
ordering of the ring $R$ must be a global ordering with
$x_i x_j > \mathrm{lm}(d_{i,j})$. The entries of the matrices $C$ and $D$ on or
below the main diagonal are ignored.

!!! note
    The conditions that assure that multiplication is associative in the
    resulting algebra are currently *not* checked.

**Examples**

```julia
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> R, (x, y) = GAlgebra(R, Singular.Matrix(R, [0 2; 0 0]),
                               Singular.Matrix(R, [0 x; 0 0]))
(Singular G-Algebra (QQ),(x,y),(dp(2),C), sgpoly{n_Q}[x, y])

julia> y*x
2*x*y + x
```

### WeylAlgebra

```julia
function WeylAlgebra(R::Union{Ring, Field}, x::Vector{String};
                     cached::Bool = true, ordering::Symbol = :degrevlex,
                     ordering2::Symbol = :comp1min, degree_bound::Int = 0)

function WeylAlgebra(R::Union{Ring, Field}, x::Matrix{String};
                     cached::Bool = true, ordering::Symbol = :degrevlex,
                     ordering2::Symbol = :comp1min, degree_bound::Int = 0)
```

Construct the ring of differential operators $\partial_1, \dots, \partial_n$
with coefficients in $R[x_1, \dots, x_n]$. In the first variant, the
differential operators are named by simply appending the letter "d" to each of
the strings in `x`. The second variant takes the names of the $x_i$ from the
first row of the matrix and the names of the $\partial_i$ from the second row
of the matrix.

**Examples**

```julia
julia> R, (x, y, dx, dy) = WeylAlgebra(ZZ, ["x", "y"])
(Singular Weyl Algebra (ZZ),(x,y,dx,dy),(dp(4),C), sweylpoly{n_Z}[x, y, dx, dy])

julia> (dx*x, dx*y, dy*x, dy*y)
(x*dx + 1, y*dx, x*dy, y*dy + 1)
```

### ExteriorAlgebra

```julia
ExteriorAlgebra(R::Union{Ring, Field}, x::Vector{String}; cached::Bool = true,
                ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min)
```

Construct the Exterior algebra as the quotient of the G-algebra with relations
$x_j x_i = -x_i x_j$ by the relations $x_i^2 = 0$.

!!! note
    We currently require at least two variables for this algebra.

**Examples**

```julia
julia> R, (x, y) = ExteriorAlgebra(ZZ, ["x", "y"])
(Singular Exterior Algebra Quotient Ring (ZZ),(x,y),(dp(2),C), sextpoly{n_Z}[x, y])

julia> x*y + y*x
0
```

### FreeAlgebra

```julia
FreeAlgebra(R::Union{Ring, Field}, x::Vector{String}, degree_bound::Int;
            cached::Bool = true, ordering = :degrevlex, ordering2::Symbol = :comp1min)
```

Construct the free associative algebra $R \langle x_1,\dots,x_n \rangle$.

!!! note
    Since this uses the LETTERPLACE backend, the `degree_bound`, which is the
    maximum length on any monomial word in the algebra, must be specified.
    Multiplication is checked and throws when the resulting degree exceeds this
    bound.

**Examples**

```julia
julia> R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 5)
(Singular letterplace Ring (QQ),(x,y,x,y,x,y,x,y,x,y),(dp(10),C,L(7)), slppoly{n_Q}[x, y])

julia> (x*y)^2
x*y*x*y

julia> (x*y)^3
    singular error: degree bound of Letterplace ring is 5, but at least 6 is needed for this multiplication
```

## Polynomial Term Iterators

For `GAlgebra`, `WeylAlgebra` and `ExteriorAlgebra`, the polynomials elements
can be and are represented using commutative data structures, and the function
`exponent_vectors` is repurposed for access to the individual exponents.

**Examples**

```julia
julia> R, (x, y, z, w) = ExteriorAlgebra(QQ, ["x", "y", "z", "w"])
(Singular Exterior Algebra Quotient Ring (QQ),(x,y,z,w),(dp(4),C), sextpoly{n_Q}[x, y, z, w])

julia> p = (1 + x + y*z)*(1 + w + x*w)
x*y*z*w + y*z*w + y*z + 2*x*w + x + w + 1

julia> show(collect(exponent_vectors(p)))
[[1, 1, 1, 1], [0, 1, 1, 1], [0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]]
```

For `FreeAlgebra`, the function `exponent_vectors` is undefined on polynomial
elements and replaced by `exponent_words` which reads off in order the indices
of the variables in a monomial. Also, the monomials for the `MPolyBuildCtx` are
specified by these exponent words. Other than these differences the term
iterators have the same behavior as in the commutative case.

**Examples**

```julia
julia> R, (x, y, z) = FreeAlgebra(QQ, ["x", "y", "z"], 6)
(Singular letterplace Ring (QQ),(x,y,z,x,y,z,x,y,z,x,y,z,x,y,z,x,y,z),(dp(18),C,L(7)), slppoly{n_Q}[x, y, z])

julia> p = (1 + x*z + y)^2
x*z*x*z + x*z*y + y*x*z + y^2 + 2*x*z + 2*y + 1

julia> show(collect(coefficients(p)))
n_Q[1, 1, 1, 1, 2, 2, 1]

julia> show(collect(monomials(p)))
slppoly{n_Q}[x*z*x*z, x*z*y, y*x*z, y*y, x*z, y, 1]

julia> show(collect(terms(p)))
slppoly{n_Q}[x*z*x*z, x*z*y, y*x*z, y*y, 2*x*z, 2*y, 1]

julia> show(collect(exponent_words(p)))
[[1, 3, 1, 3], [1, 3, 2], [2, 1, 3], [2, 2], [1, 3], [2], Int64[]]

julia> B = MPolyBuildCtx(R)
Builder for a polynomial in Singular letterplace Ring (QQ),(x,y,z,x,y,z,x,y,z,x,y,z,x,y,z,x,y,z),(dp(18),C,L(7))

julia> push_term!(B, QQ(2), [3,2,1,3]);

julia> push_term!(B, QQ(-1), Int[]);

julia> finish(B)
2*z*y*x*z - 1
```


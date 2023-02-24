```@meta
CurrentModule = Singular
```

# Noncommutative algebras

Singular.jl allows the creation of various noncommutative algebras over any of
the coefficient rings described above. The constructors of the parent objects
and elements thereof are given in the following table.

 Constructor   | Element type    | Parent type       | SINGULAR kernel subsystem
---------------|-----------------|-------------------|--------------------------
GAlgebra       | `spluralg{T}`   | `PluralRing{T}`   | PLURAL
WeylAlgebra    | `spluralg{T}`   | `PluralRing{T}`   | PLURAL
FreeAlgebra    | `slpalg{T}`     | `LPRing{T}`       | LETTERPLACE

These types are parameterized by the type of elements in the coefficient ring
of the algebra. All noncommutative algebra element types belong directly to the
abstract type `AbstractAlgebra.NCRingElem` and all the noncommuative algebra
parent object types belong to the abstract type `AbstractAlgebra.NCRing`.
The following union types cover all of Singular polynomial rings and algebras.

```julia
const PolyRingUnion{T} = Union{PolyRing{T}, PluralRing{T}, LPRing{T}} where T <: Nemo.RingElem

const SPolyUnion{T} = Union{spoly{T}, spluralg{T}, slpalg{T}} where T <: Nemo.RingElem
```

## Constructors

All constructors returns a tuple, $R, x$ consisting of a parent object $R$ and
an array $x$ of variables from which elements of the algebra can be constructed.

For constructors taking an ordering, two orderings can be specified by symbol,
one for term ordering, and a second one for ordering of module components. The
first ordering can also be specified by a non-symbol as with `polynomial_ring`,
in which case the second ordering is ignored.

By default there will only be one parent object in the system for each
combination of arguments. This is accomplished by making use of a global cache.
If this is not the desired behaviour `cached = false` may be passed.

### GAlgebra

```julia
GAlgebra(R::PolyRing{T}, C::smatrix{spoly{T}}, D::smatrix{spoly{T}};
         cached::Bool = true) where T <: Nemo.RingElem
```

Construct the G-algebra from a commutative polynomial ring $R$ and matrices $C$,
$D$ over $R$. If the variables of $R$ are $x_1,\dots,x_n$, then the noncommutative
algebra is constructed with relations $x_j x_i = c_{i,j} x_i x_j + d_{i,j}$ for
$1 \le i < j \le n$. The $c_{i,j}$ must be nonzero constant polynomials and the
relations $x_i x_j > \mathrm{lm}(d_{i,j})$ must hold in the monomial ordering
of the ring $R$.

The entries of the matrices $C$ and $D$ on or below the main diagonal are
ignored. A non-matrix argument `a` for either `C` or `D` is turned into a
matrix with all relevant entries set to `a`.

!!! note
    The conditions that assure that multiplication is associative in the
    resulting algebra are currently *not* checked. The example below
    illustrates how this condition can be checked.

**Examples**

```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> G, (x, y) = GAlgebra(R, 2, Singular.Matrix(R, [0 x; 0 0]))
(Singular G-Algebra (QQ),(x,y),(dp(2),C), spluralg{n_Q}[x, y])

julia> y*x
2*x*y + x
```

Associativity can be checked via [Interpreter Functionality](@ref).

```
julia> iszero(Singular.LibNctools.ndcond(G))
true
```

The construction of a GR-algebra proceeds by taking the quotient of a G-algebra
by a two-sided ideal. Continuing with the above example:

```
julia> I = Ideal(G, [x^2 + y^2], twosided = true)
Singular two-sided ideal over Singular G-Algebra (QQ),(x,y),(dp(2),C) with generators (x^2 + y^2)

julia> Q, (x, y) = QuotientRing(G, std(I))
(Singular G-Algebra Quotient Ring (QQ),(x,y),(dp(2),C), spluralg{n_Q}[x, y])
```

### WeylAlgebra

```julia
function WeylAlgebra(R::Union{Ring, Field}, s::Union{Vector{String}, Vector{Symbol}};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)

function WeylAlgebra(R::Union{Ring, Field}, s::Union{Matrix{String}, Matrix{Symbol}};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
```

Construct the ring of differential operators $\partial_1, \dots, \partial_n$
with coefficients in $R[x_1, \dots, x_n]$. In the first variant, the
differential operators are named by simply appending the letter "d" to each of
the strings in `x`. The second variant takes the names of the $x_i$ from the
first row of the matrix and the names of the $\partial_i$ from the second row
of the matrix. Note that the functionality of this constructor can be achieved
with the `GAlgebra` constructor: it is provided only for convenience. Note also
that due to the ordering constraint on G-algebras, the orderings `:neglex`,
`:negdeglex`, `:negdevrevlex` are excluded.

**Examples**

```julia
julia> R, (x, y, dx, dy) = WeylAlgebra(ZZ, ["x", "y"])
(Singular G-Algebra (ZZ),(x,y,dx,dy),(dp(4),C), spluralg{n_Z}[x, y, dx, dy])

julia> (dx*x, dx*y, dy*x, dy*y)
(x*dx + 1, y*dx, x*dy, y*dy + 1)
```

The ideals of G-algebras are left ideals by default.

```julia
julia> R, (x1, x2, x3, d1, d2, d3) = WeylAlgebra(QQ, ["x1" "x2" "x3"; "d1" "d2" "d3"])
(Singular G-Algebra (QQ),(x1,x2,x3,d1,d2,d3),(dp(6),C), spluralg{n_Q}[x1, x2, x3, d1, d2, d3])

julia> gens(std(Ideal(R, [x1^2*d2^2 + x2^2*d3^2, x1*d2 + x3])))
7-element Vector{spluralg{n_Q}}:
 x1*d2 + x3
 x3^2
 x2*x3 - x1
 x1*x3
 x2^2
 x1*x2
 x1^2
```

### FreeAlgebra

```julia
FreeAlgebra(R::Field, s::Union{Vector{String}, Vector{Symbol}}, degree_bound::Int;
            ordering = :degrevlex, ordering2::Symbol = :comp1min,
            cached::Bool = true)
```

Construct the free associative algebra $R \langle x_1,\dots,x_n \rangle$. The
ordering must be global.

!!! note
    Since this uses the LETTERPLACE backend, the `degree_bound`, which is the
    maximum length on any monomial word in the algebra, must be specified.
    Multiplication is checked and throws when the resulting degree exceeds this
    bound.

**Examples**

```julia
julia> R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 5)
(Singular letterplace Ring (QQ),(x,y,x,y,x,y,x,y,x,y),(dp(10),C,L(3)), slpalg{n_Q}[x, y])

julia> (x*y)^2
x*y*x*y

julia> (x*y)^3
ERROR: degree bound of Letterplace ring is 5, but at least 6 is needed for this multiplication
```

The ideals are two-sided by default for this algebra, and there is currently no
possibility of constructing one-sided ideals.

```julia
julia> R, (x, y, z) = FreeAlgebra(QQ, ["x", "y", "z"], 4)
(Singular letterplace Ring (QQ),(x,y,z,x,y,z,x,y,z,x,y,z),(dp(12),C,L(3)), slpalg{n_Q}[x, y, z])

julia> gens(std(Ideal(R, [x*y + y*z, x*x + x*y - y*x - y*y])))
8-element Vector{slpalg{n_Q}}:
 x*y + y*z
 x^2 - y*x - y^2 - y*z
 y^3 + y*z*y - y^2*z - y*z^2
 y^2*x + y*z*x + y^2*z + y*z^2
 y^2*z*y + y*z^2*y - y^2*z^2 - y*z^3
 y*z*y^2 + y*z^2*y - y*z*y*z - y*z^3
 y^2*z*x + y*z^2*x + y^2*z^2 + y*z^3
 y*z*y*x + y*z^2*x + y*z*y*z + y*z^3
```

## Term Iterators

For `GAlgebra` and `WeylAlgebra`, the elements can be and are
represented using commutative data structures, and the function
`exponent_vectors` is repurposed for access to the individual exponents.

**Examples**

```julia
julia> R, (x, y, dx, dy) = WeylAlgebra(QQ, ["x", "y"])
(Singular G-Algebra (QQ),(x,y,dx,dy),(dp(4),C), spluralg{n_Q}[x, y, dx, dy])

julia> p = (dx + dy)*(x + y)
x*dx + y*dx + x*dy + y*dy + 2

julia> show(collect(exponent_vectors(p)))
[[1, 0, 1, 0], [0, 1, 1, 0], [1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 0, 0]]
```

For `FreeAlgebra`, the function `exponent_vectors` is undefined on
elements and replaced by `exponent_words` which reads off in order the indices
of the variables in a monomial. Also, the monomials for the `MPolyBuildCtx` are
specified by these exponent words. Other than these differences the term
iterators have the same behavior as in the commutative case.

**Examples**

```julia
julia> R, (x, y, z) = FreeAlgebra(QQ, ["x", "y", "z"], 6)
(Singular letterplace Ring (QQ),(x,y,z,x,y,z,x,y,z,x,y,z,x,y,z,x,y,z),(dp(18),C,L(3)), slpalg{n_Q}[x, y, z])

julia> p = (1 + x*z + y)^2
x*z*x*z + x*z*y + y*x*z + y^2 + 2*x*z + 2*y + 1

julia> show(collect(coefficients(p)))
n_Q[1, 1, 1, 1, 2, 2, 1]

julia> show(collect(monomials(p)))
slpalg{n_Q}[x*z*x*z, x*z*y, y*x*z, y^2, x*z, y, 1]

julia> show(collect(terms(p)))
slpalg{n_Q}[x*z*x*z, x*z*y, y*x*z, y^2, 2*x*z, 2*y, 1]

julia> show(collect(exponent_words(p)))
[[1, 3, 1, 3], [1, 3, 2], [2, 1, 3], [2, 2], [1, 3], [2], Int64[]]

julia> B = MPolyBuildCtx(R)
Builder for a polynomial in Singular letterplace Ring (QQ),(x,y,z,x,y,z,x,y,z,x,y,z,x,y,z,x,y,z),(dp(18),C,L(3))

julia> push_term!(B, QQ(2), [3,2,1,3]);

julia> push_term!(B, QQ(-1), Int[]);

julia> finish(B)
2*z*y*x*z - 1
```


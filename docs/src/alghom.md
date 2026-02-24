```@meta
CurrentModule = Singular
DocTestSetup = quote
  using Singular
end
```

# Algebra Homomorphisms

Singular.jl allows the creation of algebra homomorphisms of Singular polynomial rings
over Nemo/Singular coefficient rings.

The default algebra homomorphism type in Singular.jl is the Singular `SAlgHom` type.

Additionally, a special type for the identity homomorphism has been implemented.
The type in Singular.jl for the latter is `SIdAlgHom`.

All algebra homomorphism types belong directly to the abstract type `AbstractAlgebraHomomorphism{T}`.

## Algebra Homomorphism functionality

### Constructors

Given two Singular polynomial rings $R$ and $S$ over the same base ring, the following constructors are available for creating algebra homomorphisms.

```@docs
AlgebraHomomorphism(D::PolyRing, C::PolyRing, V::Vector)
```

```@docs
IdentityAlgebraHomomorphism(D::PolyRing)
```

**Examples**

```jldoctest
julia> L = FiniteField(3, 2, "a")
(Finite field of characteristic 3 and degree 2, a)

julia> R, (x, y, z, w) = polynomial_ring(L[1], ["x", "y", "z", "w"];
                                    ordering=:negdegrevlex)
(Singular polynomial ring (9,a),(x,y,z,w),(ds(4),C), spoly{n_GF}[x, y, z, w])

julia> S, (a, b, c) = polynomial_ring(L[1], ["a", "b", "c"];
                                    ordering=:degrevlex)
(Singular polynomial ring (9,a),(a@1,b,c),(dp(3),C), spoly{n_GF}[a, b, c])

julia> V = [a, a + b^2, b - c, c + b]
4-element Vector{spoly{n_GF}}:
 a
 b^2 + a
 b + a^4*c
 b + c

julia> f = AlgebraHomomorphism(R, S, V)
Algebra homomorphism
  from Singular polynomial ring (9,a),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (9,a),(a@1,b,c),(dp(3),C)
Defining equations: spoly{n_GF}[a, b^2 + a, b + a^4*c, b + c]
```

### Operating on objects

It is possible to act on polynomials and ideals via algebra homomorphisms.

**Examples**

```jldoctest
julia> R, (x, y, z, w) = polynomial_ring(ZZ, ["x", "y", "z", "w"];
                                    ordering=:negdegrevlex)
(Singular polynomial ring (ZZ),(x,y,z,w),(ds(4),C), spoly{n_Z}[x, y, z, w])

julia> S, (a, b, c) = polynomial_ring(ZZ, ["a", "b", "c"];
                                    ordering=:degrevlex)
(Singular polynomial ring (ZZ),(a,b,c),(dp(3),C), spoly{n_Z}[a, b, c])

julia> V = [a, a + b^2, b - c, c + b]
4-element Vector{spoly{n_Z}}:
 a
 b^2 + a
 b - c
 b + c

julia> f = AlgebraHomomorphism(R, S, V)
Algebra homomorphism
  from Singular polynomial ring (ZZ),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (ZZ),(a,b,c),(dp(3),C)
Defining equations: spoly{n_Z}[a, b^2 + a, b - c, b + c]

julia> id  = IdentityAlgebraHomomorphism(S)
Identity algebra homomorphism
  from Singular polynomial ring (ZZ),(a,b,c),(dp(3),C)
  to Singular polynomial ring (ZZ),(a,b,c),(dp(3),C)
Defining equations: spoly{n_Z}[a, b, c]

julia> J = Ideal(R, [x, y^3])
Singular ideal over Singular polynomial ring (ZZ),(x,y,z,w),(ds(4),C) with generators (x, y^3)

julia> p = x + y^3 + z*w
x + z*w + y^3

julia> K = f(J)
Singular ideal over Singular polynomial ring (ZZ),(a,b,c),(dp(3),C) with generators (a, b^6 + 3*a*b^4 + 3*a^2*b^2 + a^3)

julia> q = f(p)
b^6 + 3*a*b^4 + 3*a^2*b^2 + a^3 + b^2 - c^2 + a```

### Composition

```@docs
compose(f::SAlgHom, g::SAlgHom)
```

A short command for the composition of $f$ and $g$ is `f*g`, which is the same as
`compose(f, g)`.

**Examples**

```jldoctest
julia> R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"];
                                    ordering=:negdegrevlex)
(Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C), spoly{n_Q}[x, y, z, w])

julia> S, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"];
                                    ordering=:degrevlex)
(Singular polynomial ring (QQ),(a,b,c),(dp(3),C), spoly{n_Q}[a, b, c])

julia> V = [a, a + b^2, b - c, c + b]
4-element Vector{spoly{n_Q}}:
 a
 b^2 + a
 b - c
 b + c

julia> W = [x^2, x + y + z, z*y]
3-element Vector{spoly{n_Q}}:
 x^2
 x + y + z
 y*z

julia> f = AlgebraHomomorphism(R, S, V)
Algebra homomorphism
  from Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (QQ),(a,b,c),(dp(3),C)
Defining equations: spoly{n_Q}[a, b^2 + a, b - c, b + c]

julia> g = AlgebraHomomorphism(S, R, W)
Algebra homomorphism
  from Singular polynomial ring (QQ),(a,b,c),(dp(3),C)
  to Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
Defining equations: spoly{n_Q}[x^2, x + y + z, y*z]

julia> idR  = IdentityAlgebraHomomorphism(R)
Identity algebra homomorphism
  from Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
Defining equations: spoly{n_Q}[x, y, z, w]

julia> h1 = f*g
Algebra homomorphism
  from Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
Defining equations: spoly[x^2, 2*x^2 + 2*x*y + y^2 + 2*x*z + 2*y*z + z^2, x + y + z - y*z, x + y + z + y*z]

julia> h2 = idR*f
Algebra homomorphism
  from Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (QQ),(a,b,c),(dp(3),C)
Defining equations: spoly{n_Q}[a, b^2 + a, b - c, b + c]

julia> h3 = g*idR
Algebra homomorphism
  from Singular polynomial ring (QQ),(a,b,c),(dp(3),C)
  to Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
Defining equations: spoly{n_Q}[x^2, x + y + z, y*z]

julia> h4 = idR*idR
Identity algebra homomorphism
  from Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
Defining equations: spoly{n_Q}[x, y, z, w]
```

### Preimages

```@docs
preimage(f::SAlgHom, I::sideal)
```

```@docs
preimage(f::SIdAlgHom, I::sideal)
```

```@docs
kernel(f::SIdAlgHom)
```

```@docs
kernel(f::SAlgHom)
```

**Examples**

```jldoctest
julia> R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"];
                                    ordering=:negdegrevlex)
(Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C), spoly{n_Q}[x, y, z, w])

julia> S, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"];
                                    ordering=:degrevlex)
(Singular polynomial ring (QQ),(a,b,c),(dp(3),C), spoly{n_Q}[a, b, c])

julia> I = Ideal(S, [a, a + b^2, b - c, c + b])
Singular ideal over Singular polynomial ring (QQ),(a,b,c),(dp(3),C) with generators (a, b^2 + a, b - c, b + c)

julia> f = AlgebraHomomorphism(R, S, gens(I))
Algebra homomorphism
  from Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C)
  to Singular polynomial ring (QQ),(a,b,c),(dp(3),C)
Defining equations: spoly{n_Q}[a, b^2 + a, b - c, b + c]

julia> idS  = IdentityAlgebraHomomorphism(S)
Identity algebra homomorphism
  from Singular polynomial ring (QQ),(a,b,c),(dp(3),C)
  to Singular polynomial ring (QQ),(a,b,c),(dp(3),C)
Defining equations: spoly{n_Q}[a, b, c]

julia> P1 = preimage(f, I)
Singular ideal over Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C) with generators (x, y, z, w)

julia> P2 = preimage(idS, I)
Singular ideal over Singular polynomial ring (QQ),(a,b,c),(dp(3),C) with generators (a, b^2 + a, b - c, b + c)

julia> K1 = kernel(f)
Singular ideal over Singular polynomial ring (QQ),(x,y,z,w),(ds(4),C) with generators (4*x - 4*y + z^2 + 2*z*w + w^2)

julia> K2 = kernel(idS)
Singular ideal over Singular polynomial ring (QQ),(a,b,c),(dp(3),C) with generators (0)
```


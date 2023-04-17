```@meta
CurrentModule = Singular
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

```julia
L = FiniteField(3, 2, "a")

R, (x, y, z, w) = polynomial_ring(L[1], ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)

S, (a, b, c) = polynomial_ring(L[1], ["a", "b", "c"];
                             ordering=:degrevlex)

V = [a, a + b^2, b - c, c + b]

f = AlgebraHomomorphism(R, S, V)
```

### Operating on objects

It is possible to act on polynomials and ideals via algebra homomorphisms.

**Examples**

```
R, (x, y, z, w) = polynomial_ring(Nemo.ZZ, ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)

S, (a, b, c) = polynomial_ring(Nemo.ZZ, ["a", "b", "c"];
                             ordering=:degrevlex)

V = [a, a + b^2, b - c, c + b]

f = AlgebraHomomorphism(R, S, V)

id  = IdentityAlgebraHomomorphism(S)


J = Ideal(R, [x, y^3])

p = x + y^3 + z*w

K = f(J)

q = f(p)
```

### Composition

```@docs
compose(f::SAlgHom, g::SAlgHom)
```

A short command for the composition of $f$ and $g$ is `f*g`, which is the same as
`compose(f, g)`.

**Examples**

```
R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)

S, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"];
                             ordering=:degrevlex)

V = [a, a + b^2, b - c, c + b]

W = [x^2, x + y + z, z*y]

f = AlgebraHomomorphism(R, S, V)

g = AlgebraHomomorphism(S, R, W)

idR  = IdentityAlgebraHomomorphism(R)

h1 = f*g

h2 = idR*f

h3 = g*idR

h4 = idR*idR
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

```
R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)

S, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"];
                             ordering=:degrevlex)

I = Ideal(S, [a, a + b^2, b - c, c + b])

f = SAlgebraHomomorphism(R, S, gens(I))

idS  = IdentityAlgebraHomomorphism(S)

P1 = preimage(f, I)

P2 = preimage(idS, I)

K1 = kernel(f)

K2 = preimage(idS)
```


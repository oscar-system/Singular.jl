```@meta
CurrentModule = Singular
```

# Algebra Homomorphisms

Singular.jl allows the creation of algebra homomorphisms of Singular polynomial rings over the fields `QQ`, `Fp` and finite extensions of `Fp`. 

The default algebra homomorphism type in Singular.jl is the Singular `SAlgHom` type.

Additionally, a special type for the identity homomorphism has been implemented. 
The type in Singular.jl for the latter is `SIdAlgHom`. 

All algebra homomorphism types belong directly to the abstract type `AbstractAlgebraHomomorphism{T}`.

## Algebra Homomorphism functionality 

### Constructors

Given two Singular polynomial rings $R$ and $S$ over the same base field $K$, the following constructors are available for creating ideals.

```@docs
AlgebraHomomorphism(D::PolyRing, C::PolyRing, I::sideal)
```

```@docs
IdentityAlgebraHomomorphism(D::PolyRing)
```

**Examples**

```julia
L = FiniteField(3, 2, String("a"))

R, (x, y, z, w) = PolynomialRing(L[1], ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)

S, (a, b, c) = PolynomialRing(L[1], ["a", "b", "c"];
                             ordering=:degrevlex)

I = Ideal(S, [a, a + b^2, b - c, c + b])

f = AlgebraHomomorphism(R, S, I)
```

### Operating on objects

It is possible to act on polynomials and ideals via algebra homomorphisms. 

**Examples**

```
R, (x, y, z, w) = PolynomialRing(Singular.Fp(3), ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)
   
S, (a, b, c) = PolynomialRing(Singular.Fp(3), ["a", "b", "c"];
                             ordering=:degrevlex)
   
I = Ideal(S, [a, a + b^2, b - c, c + b])
   
f = AlgebraHomomorphism(R, S, I)
   
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
R, (x, y, z, w) = PolynomialRing(Singular.QQ, ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)

S, (a, b, c) = PolynomialRing(Singular.QQ, ["a", "b", "c"];
                             ordering=:degrevlex)

I = Ideal(S, [a, a + b^2, b - c, c + b])

K = Ideal(R, [x^2, x + y + z, z*y])

L = Ideal(R, [x^2, 2*x^2 + 2*x*y + y^2 + 2*x*z + 2*y*z + z^2, x + y + z - y*z,
                x + y + z + y*z])

f = AlgebraHomomorphism(R, S, I)

g = AlgebraHomomorphism(S, R, K)

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
kernel(f::SAlgHom)
```

**Examples**

```
R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)
   
S, (a, b, c) = PolynomialRing(QQ, ["a", "b", "c"];
                             ordering=:degrevlex)
   
I = Ideal(S, [a, a + b^2, b - c, c + b])
   
f = SAlgebraHomomorphism(R, S, I)
   
idS  = IdentityAlgebraHomomorphism(S)

P = preimage(f, I)

K = kernel(f)
```


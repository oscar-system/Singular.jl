```@meta
CurrentModule = Singular
```

# Quotient Rings

Quotient rings $Q = R/I$ in Singular.jl are constructed with the constructor
`QuotientRing(R, I)`. The input ideal $I$ to the constructor must be a Groebner
basis. The $R$-ideal $I$ may be recovered as `quotient_ideal(Q)`.

```@docs
is_quotient_ring(R::PolyRingUnion)
quotient_ideal(Q::PolyRing{T}) where T <: Nemo.RingElem
```

**Examples**

```julia
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> is_quotient_ring(R)
false

julia> Q1, (x, y) = QuotientRing(R, Ideal(R, x^2+y^2));

julia> is_quotient_ring(Q1)
true

julia> quotient_ideal(Q1)
Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x^2 + y^2)

julia> Q2, (x, y) = QuotientRing(Q1, std(Ideal(Q1, x*y)));

julia> quotient_ideal(Q2)
Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x*y, y^3, x^2 + y^2)

julia> base_ring(quotient_ideal(Q1)) == base_ring(quotient_ideal(Q2))
true
```

```@meta
CurrentModule = Singular
```

# Interpreter functionality

## Library procedures

## Global variables

The global variables `degBound` and `multBound` can be used in a local fashion.
As with any global variable, their usage should be accompanied with caution.

```@docs
with_degBound(f, degb::Integer)
```

```@docs
with_multBound(f, degb::Integer)
```
**Examples**

```julia
julia> r, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=ordering_ds());

julia> i = Ideal(r, [x^7+y^7+z^6,x^6+y^8+z^7,x^7+y^5+z^8,x^2*y^3+y^2*z^3+x^3*z^2,x^3*y^2+y^3*z^2+x^2*z^3]);

julia> degree(std(i))   # default behaviour of no multiplicity bound
(0, 86)

julia> with_multBound(100) do
           # run with a multiplicity bound of 100
           return degree(std(i))
       end
(0, 98)

julia> degree(std(i))   # back to default behaviour
(0, 86)

julia> gens(std(i))
11-element Vector{spoly{n_Q}}:
 x^3*y^2 + y^3*z^2 + x^2*z^3
 x^2*y^3 + x^3*z^2 + y^2*z^3
 y^5 + x^7 + z^8
 x^6 + z^7 + y^8
 x^4*z^2 - y^4*z^2 - x^2*y*z^3 + x*y^2*z^3
 z^6 + x^7 + y^7
 y^4*z^3 - y^3*z^4 - x^2*z^5 - x^9
 x^3*y*z^4 - x^2*y^2*z^4 + x*y^3*z^4 - y^4*z^4 + x^3*z^5 - x^2*y*z^5
 x^3*z^5
 x^2*y*z^5 + y^3*z^5 + x^2*z^6
 x*y^3*z^5

julia> gens(with_degBound(5) do; return std(i); end)
5-element Vector{spoly{n_Q}}:
 x^3*y^2 + y^3*z^2 + x^2*z^3
 x^2*y^3 + x^3*z^2 + y^2*z^3
 y^5 + x^7 + z^8
 x^6 + z^7 + y^8
 z^6 + x^7 + y^7
```


"""
    Singular.randseed!([seed::Integer])

Reseed Singular's global RNG with `seed`.

The given `seed` must be a non-zero integer.
When `seed` is not specified, a random seed is generated from Julia's global RNG.
"""
randseed!(seed::Integer) = randseed!(seed % Cint)
randseed!() = randseed!(rand(UInt128))
function randseed!(seed::Cint)
    Singular.libSingular.set_randomseed(seed)
    nothing
end

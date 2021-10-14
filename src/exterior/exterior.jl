export pexterior,
       ExteriorAlgebra

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::ExteriorAlgebra{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::pexterior) = base_ring(parent(p))

elem_type(::Type{ExteriorAlgebra{T}}) where T <: Nemo.RingElem = pexterior{T}

parent_type(::Type{pexterior{T}}) where T <: Nemo.RingElem = ExteriorAlgebra{T}

@doc Markdown.doc"""
    nvars(R::ExteriorAlgebra)
> Return the number of variables in the given Exterior algebra.
"""
nvars(R::ExteriorAlgebra) = Int(libSingular.rVar(R.ptr))

@doc Markdown.doc"""
    has_global_ordering(R::ExteriorAlgebra)
> Return `true` if the given algebra has a global ordering, i.e. if $1 < x$ for
> each variable $x$ in the ring. This include `:lex`, `:deglex` and `:degrevlex`
> orderings.
"""
has_global_ordering(R::ExteriorAlgebra) = Bool(libSingular.rHasGlobalOrdering(R.ptr))

@doc Markdown.doc"""
    has_mixed_ordering(R::ExteriorAlgebra)
> Return `true` if the given algebra has a mixed ordering, i.e. if $1 < x_i$ for
> a variable $x_i$ and $1>x_j$ for another variable $x_j$.
"""
has_mixed_ordering(R::ExteriorAlgebra) = Bool(libSingular.rHasMixedOrdering(R.ptr))

@doc Markdown.doc"""
    has_local_ordering(R::ExteriorAlgebra)
> Return `true` if the given algebra has a local ordering, i.e. if $1 > x$ for
> all variables $x$.
"""
function has_local_ordering(R::ExteriorAlgebra)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

@doc Markdown.doc"""
    isquotient_ring(R::ExteriorAlgebra)
> Return `true` if the given algebra is the quotient of a polynomial ring with
> a non - zero ideal.
"""
isquotient_ring(R::ExteriorAlgebra) = Bool(Singular.libSingular.rIsQuotientRing(R.ptr))

@doc Markdown.doc"""
    characteristic(R::ExteriorAlgebra)
> Return the characteristic of the Exterior algebra, i.e. the characteristic of the
> coefficient ring.
"""
characteristic(R::ExteriorAlgebra) = Int(libSingular.rChar(R.ptr))

function gens(R::ExteriorAlgebra)
   n = nvars(R)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::ExteriorAlgebra, i::Int)
   return R(libSingular.rGetVar(Cint(i), R.ptr))
end

@doc Markdown.doc"""
    symbols(R::ExteriorAlgebra)
> Return symbols for the generators of the Exterior algebra $R$.
"""
function symbols(R::ExteriorAlgebra)
   io = IOBuffer();
   symbols = Array{Symbol,1}(undef, 0)
   for g in gens(R)
      show(io,g)
      push!(symbols, Symbol(String(take!(io))))
   end
   return symbols
end

ordering(R::ExteriorAlgebra) = R.ord

@doc Markdown.doc"""
    degree_bound(R::ExteriorAlgebra)
> Return the internal degree bound in each variable, enforced by Singular. This is the
> largest positive value any degree can have before an overflow will occur. This
> internal bound may be higher than the bound requested by the user via the
> `degree_bound` parameter of the `ExteriorAlgebra` constructor.
"""
function degree_bound(R::ExteriorAlgebra)
   return Int(libSingular.rBitmask(R.ptr))
end

zero(R::ExteriorAlgebra) = R()

one(R::ExteriorAlgebra) = R(1)


function Base.hash(p::pexterior{T}, h::UInt) where T <: Nemo.RingElem
   v = 0xaf708b07f940b4d2%UInt
   v = xor(hash(collect(exponent_vectors(p)), h), v)
   for c in coefficients(p)
      v = xor(hash(c, h), v)
      v = (v << 1) | (v >> (sizeof(Int)*8 - 1))
   end
   return v
end

###############################################################################
#
#   String I/O
#
###############################################################################          
function show(io::IO, R::ExteriorAlgebra)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
     print(io, "Singular Exterior Algebra Quotient Ring ", s)
   else
     print(io, "Singular Exterior Algebra ", s)
   end
end

show_minus_one(::Type{pexterior{T}}) where T <: Nemo.RingElem = show_minus_one(T)

needs_parentheses(x::pexterior) = length(x) > 1

isnegative(x::pexterior) = isconstant(x) && !iszero(x) && isnegative(coeff(x, 0))

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{pexterior{T}}, ::Type{pexterior{T}}) where T <: Nemo.RingElem = pexterior{T}

function promote_rule(::Type{pexterior{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? pexterior{T} : Union{}
end

###############################################################################
#
#   Build context
#
###############################################################################

function MPolyBuildCtx(R::ExteriorAlgebra)
   return MPolyBuildCtx(R, R().ptr)
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::ExteriorAlgebra)()
   T = elem_type(base_ring(R))
   z = pexterior{T}(R)
   z.parent = R
   return z
end

function (R::ExteriorAlgebra)(n::Int)
   T = elem_type(base_ring(R))
   z = pexterior{T}(R, n)
   z.parent = R
   return z
end

function (R::ExteriorAlgebra)(n::Integer)
   T = elem_type(base_ring(R))
   z = pexterior{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::ExteriorAlgebra)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = pexterior{T}(R, ptr)
   z.parent = R
   return z
end

function (R::ExteriorAlgebra)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   z = pexterior{T}(R, n)
   z.parent = R
   return z
end

function (R::ExteriorAlgebra{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Exterior algebra")
   z = pexterior{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::ExteriorAlgebra{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::ExteriorAlgebra)(p::pexterior)
   parent(p) != R && error("Unable to coerce")
   return p
end

function(R::ExteriorAlgebra)(n::libSingular.number_ptr)
    return R.base_ring(n)
end
###############################################################################
#
#   ExteriorAlgebra constructor
#
###############################################################################

function ExteriorAlgebra(R::Union{Ring, Field}, s::Array{String, 1};
      cached::Bool = true, ordering::Symbol = :degrevlex,
      ordering2::Symbol = :comp1min, degree_bound::Int = 0)
   U = [Symbol(v) for v in s]
   T = elem_type(R)
   parent_obj = ExteriorAlgebra{T}(R, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

function ExteriorAlgebra(R::Nemo.Ring, s::Array{String, 1}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)
   S = CoefficientRing(R)
   U = [Symbol(v) for v in s]
   T = elem_type(S)
   parent_obj = ExteriorAlgebra{T}(S, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

macro ExteriorAlgebra(R, s, n, o)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = ExteriorAlgebra($R, $v0; ordering=$o))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

macro ExteriorAlgebra(R, s, n)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = ExteriorAlgebra($R, $v0))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end
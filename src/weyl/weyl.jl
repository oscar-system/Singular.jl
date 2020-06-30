###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::WeylAlgebra{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::pweyl) = base_ring(parent(p))

elem_type(::Type{WeylAlgebra{T}}) where T <: Nemo.RingElem = pweyl{T}

parent_type(::Type{pweyl{T}}) where T <: Nemo.RingElem = WeylAlgebra{T}

@doc Markdown.doc"""
    nvars(R::WeylAlgebra)
> Return the number of variables in the given Weyl algebra.
"""
nvars(R::WeylAlgebra) = Int(libSingular.rVar(R.ptr))

@doc Markdown.doc"""
    has_global_ordering(R::WeylAlgebra)
> Return `true` if the given algebra has a global ordering, i.e. if $1 < x$ for
> each variable $x$ in the ring. This include `:lex`, `:deglex` and `:degrevlex`
> orderings.
"""
has_global_ordering(R::WeylAlgebra) = Bool(libSingular.rHasGlobalOrdering(R.ptr))

@doc Markdown.doc"""
    has_mixed_ordering(R::WeylAlgebra)
> Return `true` if the given algebra has a mixed ordering, i.e. if $1 < x_i$ for
> a variable $x_i$ and $1>x_j$ for another variable $x_j$.
"""
has_mixed_ordering(R::WeylAlgebra) = Bool(libSingular.rHasMixedOrdering(R.ptr))

@doc Markdown.doc"""
    has_local_ordering(R::WeylAlgebra)
> Return `true` if the given algebra has a local ordering, i.e. if $1 > x$ for
> all variables $x$.
"""
function has_local_ordering(R::WeylAlgebra)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

@doc Markdown.doc"""
    isquotient_ring(R::WeylAlgebra)
> Return `true` if the given algebra is the quotient of a polynomial ring with
> a non - zero ideal.
"""
isquotient_ring(R::WeylAlgebra) = Bool(Singular.libSingular.rIsQuotientRing(R.ptr))

@doc Markdown.doc"""
    characteristic(R::WeylAlgebra)
> Return the characteristic of the Weyl algebra, i.e. the characteristic of the
> coefficient ring.
"""
characteristic(R::WeylAlgebra) = Int(libSingular.rChar(R.ptr))

function gens(R::WeylAlgebra)
   n = nvars(R)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::WeylAlgebra, i::Int)
   return R(libSingular.rGetVar(Cint(i), R.ptr))
end

@doc Markdown.doc"""
    symbols(R::WeylAlgebra)
> Return symbols for the generators of the Weyl algebra $R$.
"""
function symbols(R::WeylAlgebra)
   io = IOBuffer();
   symbols = Array{Symbol,1}(undef, 0)
   for g in gens(R)
      show(io,g)
      push!(symbols, Symbol(String(take!(io))))
   end
   return symbols
end

ordering(R::WeylAlgebra) = R.ord

@doc Markdown.doc"""
    degree_bound(R::WeylAlgebra)
> Return the internal degree bound in each variable, enforced by Singular. This is the
> largest positive value any degree can have before an overflow will occur. This
> internal bound may be higher than the bound requested by the user via the
> `degree_bound` parameter of the `WeylAlgebra` constructor.
"""
function degree_bound(R::WeylAlgebra)
   return Int(libSingular.rBitmask(R.ptr))
end

zero(R::WeylAlgebra) = R()

one(R::WeylAlgebra) = R(1)

function Base.hash(p::pweyl{T}, h::UInt) where T <: Nemo.RingElem
   v = 0xa4c406868a7c20c5%UInt
   v = xor(hash(collect(exponent_vectors(p)), h), v)
   for c in coeffs(p)
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
function show(io::IO, R::WeylAlgebra)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
     print(io, "Singular Weyl Algebra Quotient Ring ", s)
   else
     print(io, "Singular Weyl Algebra ", s)
   end
end

show_minus_one(::Type{pweyl{T}}) where T <: Nemo.RingElem = show_minus_one(T)

needs_parentheses(x::pweyl) = length(x) > 1

isnegative(x::pweyl) = isconstant(x) && !iszero(x) && isnegative(coeff(x, 0))

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{pweyl{T}}, ::Type{pweyl{T}}) where T <: Nemo.RingElem = pweyl{T}

function promote_rule(::Type{pweyl{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? pweyl{T} : Union{}
end

###############################################################################
#
#   Build context
#
###############################################################################

function MPolyBuildCtx(R::WeylAlgebra)
   return MPolyBuildCtx(R, R().ptr)
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::WeylAlgebra)()
   T = elem_type(base_ring(R))
   z = pweyl{T}(R)
   z.parent = R
   return z
end

function (R::WeylAlgebra)(n::Int)
   T = elem_type(base_ring(R))
   z = pweyl{T}(R, n)
   z.parent = R
   return z
end

function (R::WeylAlgebra)(n::Integer)
   T = elem_type(base_ring(R))
   z = pweyl{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::WeylAlgebra)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = pweyl{T}(R, ptr)
   z.parent = R
   return z
end

function (R::WeylAlgebra)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   z = pweyl{T}(R, n)
   z.parent = R
   return z
end

function (R::WeylAlgebra{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Weyl algebra")
   z = pweyl{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::WeylAlgebra{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::WeylAlgebra)(p::pweyl)
   parent(p) != R && error("Unable to coerce")
   return p
end

function(R::WeylAlgebra)(n::libSingular.number_ptr)
    return R.base_ring(n)
end
###############################################################################
#
#   WeylAlgebra constructor
#
###############################################################################

function WeylAlgebra(R::Union{Ring, Field}, s1::Array{String, 1};
      cached::Bool = true, ordering::Symbol = :degrevlex,
      ordering2::Symbol = :comp1min, degree_bound::Int = 0)
   s = vcat(s1, ["d"*sym for sym in s1])
   U = [Symbol(v) for v in s]
   T = elem_type(R)
   parent_obj = WeylAlgebra{T}(R, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

function WeylAlgebra(R::Union{Ring, Field}, s2::Array{String, 2};
      cached::Bool = true, ordering::Symbol = :degrevlex,
      ordering2::Symbol = :comp1min, degree_bound::Int = 0)
   size(s2)[1] != 2 && error("Array of symbols should have two rows")
   s = vcat(s2[1, :], s2[2, :])
   U = [Symbol(v) for v in s]
   T = elem_type(R)
   parent_obj = WeylAlgebra{T}(R, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

function WeylAlgebra(R::Nemo.Ring, s1::Array{String, 1}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)
   s = vcat(s1, ["d"*sym for sym in s1])
   S = CoefficientRing(R)
   U = [Symbol(v) for v in s]
   T = elem_type(S)
   parent_obj = WeylAlgebra{T}(S, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

function WeylAlgebra(R::Nemo.Ring, s2::Array{String, 2}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)
   size(s2)[1] != 2 && error("Array of symbols should have two rows")
   s = vcat(s2[1, :], s2[2, :])
   S = CoefficientRing(R)
   U = [Symbol(v) for v in s]
   T = elem_type(S)
   parent_obj = WeylAlgebra{T}(S, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

macro WeylAlgebra(R, s, n, o)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = WeylAlgebra($R, $v0; ordering=$o))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

macro WeylAlgebra(R, s, n)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = WeylAlgebra($R, $v0))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

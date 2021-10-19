export sweylpoly, WeylPolyRing, WeylAlgebra

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::WeylPolyRing{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::sweylpoly) = base_ring(parent(p))

elem_type(::Type{WeylPolyRing{T}}) where T <: Nemo.RingElem = sweylpoly{T}

parent_type(::Type{sweylpoly{T}}) where T <: Nemo.RingElem = WeylPolyRing{T}

nvars(R::WeylPolyRing) = Int(libSingular.rVar(R.ptr))

has_global_ordering(R::WeylPolyRing) = Bool(libSingular.rHasGlobalOrdering(R.ptr))

has_mixed_ordering(R::WeylPolyRing) = Bool(libSingular.rHasMixedOrdering(R.ptr))

function has_local_ordering(R::WeylPolyRing)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

isquotient_ring(R::WeylPolyRing) = Bool(Singular.libSingular.rIsQuotientRing(R.ptr))

characteristic(R::WeylPolyRing) = Int(libSingular.rChar(R.ptr))

function gens(R::WeylPolyRing)
   n = nvars(R)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::WeylPolyRing, i::Int)
   return R(libSingular.rGetVar(Cint(i), R.ptr))
end

ordering(R::WeylPolyRing) = R.ord

@doc Markdown.doc"""
    degree_bound(R::WeylPolyRing)

Return the internal degree bound in each variable, enforced by Singular. This is the
largest positive value any degree can have before an overflow will occur. This
internal bound may be higher than the bound requested by the user via the
`degree_bound` parameter of the `WeylAlgebra` constructor.
"""
function degree_bound(R::WeylPolyRing)
   return Int(libSingular.rBitmask(R.ptr))
end

zero(R::WeylPolyRing) = R()

one(R::WeylPolyRing) = R(1)

function Base.hash(p::sweylpoly{T}, h::UInt) where T <: Nemo.RingElem
   v = 0xa4c406868a7c20c5%UInt
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
function show(io::IO, R::WeylPolyRing)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
     print(io, "Singular Weyl Algebra Quotient Ring ", s)
   else
     print(io, "Singular Weyl Algebra ", s)
   end
end

show_minus_one(::Type{sweylpoly{T}}) where T <: Nemo.RingElem = show_minus_one(T)

needs_parentheses(x::sweylpoly) = length(x) > 1

isnegative(x::sweylpoly) = isconstant(x) && !iszero(x) && isnegative(coeff(x, 0))

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{sweylpoly{T}}, ::Type{sweylpoly{T}}) where T <: Nemo.RingElem = sweylpoly{T}

function promote_rule(::Type{sweylpoly{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? sweylpoly{T} : Union{}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::WeylPolyRing)()
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R)
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::Int)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, ptr)
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::WeylPolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Weyl algebra")
   z = sweylpoly{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::WeylPolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::WeylPolyRing)(p::sweylpoly)
   parent(p) != R && error("Unable to coerce")
   return p
end

function(R::WeylPolyRing)(n::libSingular.number_ptr)
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
   parent_obj = WeylPolyRing{T}(R, U, ordering, cached, sym2ringorder[ordering],
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
   parent_obj = WeylPolyRing{T}(R, U, ordering, cached, sym2ringorder[ordering],
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
   parent_obj = WeylPolyRing{T}(S, U, ordering, cached, sym2ringorder[ordering],
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
   parent_obj = WeylPolyRing{T}(S, U, ordering, cached, sym2ringorder[ordering],
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

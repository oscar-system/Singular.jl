export sextpoly, ExtPolyRing, ExteriorAlgebra

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::ExtPolyRing{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::sextpoly) = base_ring(parent(p))

elem_type(::Type{ExtPolyRing{T}}) where T <: Nemo.RingElem = sextpoly{T}

parent_type(::Type{sextpoly{T}}) where T <: Nemo.RingElem = ExtPolyRing{T}

function degree_bound(R::ExtPolyRing)
   return Int(libSingular.rBitmask(R.ptr))
end

###############################################################################
#
#   String I/O
#
###############################################################################          

function show(io::IO, R::ExtPolyRing)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
     print(io, "Singular Exterior Algebra Quotient Ring ", s)
   else
     print(io, "Singular Exterior Algebra ", s)
   end
end

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{sextpoly{T}}, ::Type{sextpoly{T}}) where T <: Nemo.RingElem = sextpoly{T}

function promote_rule(::Type{sextpoly{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? sextpoly{T} : Union{}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::ExtPolyRing)()
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R)
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::Int)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, ptr)
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::ExtPolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Exterior algebra")
   z = sextpoly{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::ExtPolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::ExtPolyRing)(p::sextpoly)
   parent(p) != R && error("Unable to coerce")
   return p
end

function (R::ExtPolyRing)(n::libSingular.number_ptr)
    return R.base_ring(n)
end

###############################################################################
#
#   ExtPolyRing constructor
#
###############################################################################

function _ExteriorAlgebra(R, s::Vector{String}, ordering, ordering2, cached)
   sord = get_fancy_ordering(ordering, ordering2)
   z = ExtPolyRing{elem_type(R)}(R, Symbol.(s), sord, cached)
   return (z, gens(z))
end

function ExteriorAlgebra(R::Union{Ring, Field}, s::Vector{String};
                         ordering = :degrevlex, ordering2::Symbol = :comp1min,
                         cached::Bool = true)
   return _ExteriorAlgebra(R, s, ordering, ordering2, cached)
end

function ExteriorAlgebra(R::Nemo.Ring, s::Vector{String};
                         ordering = :degrevlex, ordering2::Symbol = :comp1min,
                         cached::Bool = true)
   R = CoefficientRing(R)
   return _ExteriorAlgebra(R, s, ordering, ordering2, cached)
end

macro ExteriorAlgebra(R, s, n, o)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = ExtPolyRing($R, $v0; ordering=$o))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

macro ExteriorAlgebra(R, s, n)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = ExtPolyRing($R, $v0))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

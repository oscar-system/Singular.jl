export sgpoly, GPolyRing, GAlgebra

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::GPolyRing{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::sgpoly) = base_ring(parent(p))

elem_type(::Type{GPolyRing{T}}) where T <: Nemo.RingElem = sgpoly{T}

parent_type(::Type{sgpoly{T}}) where T <: Nemo.RingElem = GPolyRing{T}

@doc Markdown.doc"""
    degree_bound(R::GPolyRing)

Return the internal degree bound in each variable, enforced by Singular. This is the
largest positive value any degree can have before an overflow will occur.
"""
function degree_bound(R::GPolyRing)
   return Int(libSingular.rBitmask(R.ptr))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::GPolyRing)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
     print(io, "Singular G-Algebra Quotient Ring ", s)
   else
     print(io, "Singular G-Algebra ", s)
   end
end

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{sgpoly{T}}, ::Type{sgpoly{T}}) where T <: Nemo.RingElem = sgpoly{T}

function promote_rule(::Type{sgpoly{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? sgpoly{T} : Union{}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::GPolyRing)()
   T = elem_type(base_ring(R))
   z = sgpoly{T}(R)
   z.parent = R
   return z
end

function (R::GPolyRing)(n::Int)
   T = elem_type(base_ring(R))
   z = sgpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::GPolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   z = sgpoly{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::GPolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = sgpoly{T}(R, ptr)
   z.parent = R
   return z
end

function (R::GPolyRing)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   z = sgpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::GPolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Exterior algebra")
   z = sgpoly{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::GPolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::GPolyRing)(p::sgpoly)
   parent(p) != R && error("Unable to coerce")
   return p
end

function (R::GPolyRing)(n::libSingular.number_ptr)
    return R.base_ring(n)
end

###############################################################################
#
#   GAlgebra constructor
#
###############################################################################

function GAlgebra(R::PolyRing{T}, C::smatrix{spoly{T}}, D::smatrix{spoly{T}};
                  cached::Bool = true) where T <: Nemo.RingElem
   parent_obj = GPolyRing{T}(R, C, D, R.S)
   return (parent_obj, gens(parent_obj))
end


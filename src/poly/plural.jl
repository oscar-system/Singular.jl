export spluralg, PluralRing, GAlgebra

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::PluralRing{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::spluralg) = base_ring(parent(p))

elem_type(::Type{PluralRing{T}}) where T <: Nemo.RingElem = spluralg{T}

parent_type(::Type{spluralg{T}}) where T <: Nemo.RingElem = PluralRing{T}

@doc Markdown.doc"""
    degree_bound(R::PluralRing)

Return the internal degree bound in each variable, enforced by Singular. This is the
largest positive value any degree can have before an overflow will occur.
"""
function degree_bound(R::PluralRing)
   return Int(libSingular.rBitmask(R.ptr))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::PluralRing)
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

promote_rule(::Type{spluralg{T}}, ::Type{spluralg{T}}) where T <: Nemo.RingElem = spluralg{T}

function promote_rule(::Type{spluralg{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? spluralg{T} : Union{}
end

function promote_rule(::Type{spluralg{T}}, ::Type{T}) where {T <: Nemo.RingElem}
   return spluralg{T}
end

function promote_rule(::Type{spluralg{n_RingElem{RingElemWrapper{S, T}}}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem, S}
   return spluralg{n_RingElem{RingElemWrapper{S, T}}}
end

function promote_rule(::Type{spluralg{n_FieldElem{FieldElemWrapper{S, T}}}}, ::Type{U}) where {T <: Nemo.FieldElem, U <: Nemo.FieldElem, S}
   return spluralg{n_FieldElem{FieldElemWrapper{S, T}}}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::PluralRing{T})() where T <: Nemo.RingElem
   return spluralg{T}(R)
end

function (R::PluralRing{T})(n::Int) where T <: Nemo.RingElem
   return spluralg{T}(R, n)
end

function (R::PluralRing{T})(n::Integer) where T <: Nemo.RingElem
   return spluralg{T}(R, BigInt(n))
end

function (R::PluralRing{T})(n::n_Z) where T <: Nemo.RingElem
   return spluralg{T}(R, base_ring(R)(n))
end

function (R::PluralRing{T})(n::Rational) where T <: Nemo.RingElem
   return spluralg{T}(R, base_ring(R)(n))
end

# take ownership of the pointer - not for general users
function (R::PluralRing{T})(ptr::libSingular.poly_ptr) where T <: Nemo.RingElem
   return spluralg{T}(R, ptr)
end

function (R::PluralRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Exterior algebra")
   return spluralg{T}(R, n)
end

function (R::PluralRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::PluralRing)(p::spluralg)
   parent(p) != R && error("Unable to coerce")
   return p
end

###############################################################################
#
#   GAlgebra constructor
#
###############################################################################

function _to_matrix(R::PolyRing, a::smatrix)
   R == base_ring(a) || error("matrix has the wrong base ring")
   return a
end

function _to_matrix(R::PolyRing, a)
   n = nvars(R)
   m = zero_matrix(R, n, n)
   for i in 1:n, j in i+1:n
      m[i,j] = R(a)
   end
   return m
end

function GAlgebra(R::PolyRing{T}, C, D; cached::Bool = true) where T <: Nemo.RingElem
   parent_obj = PluralRing{T}(R, _to_matrix(R, C), _to_matrix(R, D), R.S)
   return (parent_obj, gens(parent_obj))
end


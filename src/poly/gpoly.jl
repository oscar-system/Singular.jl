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

nvars(R::GPolyRing) = Int(libSingular.rVar(R.ptr))

has_global_ordering(R::GPolyRing) = Bool(libSingular.rHasGlobalOrdering(R.ptr))

has_mixed_ordering(R::GPolyRing) = Bool(libSingular.rHasMixedOrdering(R.ptr))

function has_local_ordering(R::GPolyRing)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

isquotient_ring(R::GPolyRing) = Bool(Singular.libSingular.rIsQuotientRing(R.ptr))

characteristic(R::GPolyRing) = Int(libSingular.rChar(R.ptr))

function gens(R::GPolyRing)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i in 1:nvars(R)]
end

function gen(R::GPolyRing, i::Int)
   return R(libSingular.rGetVar(Cint(i), R.ptr))
end

function symbols(R::GPolyRing)
   return R.S
end

ordering(R::GPolyRing) = R.ord

function degree_bound(R::GPolyRing)
   return Int(libSingular.rBitmask(R.ptr))
end

zero(R::GPolyRing) = R()

one(R::GPolyRing) = R(1)

function Base.hash(p::sgpoly{T}, h::UInt) where T <: Nemo.RingElem
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


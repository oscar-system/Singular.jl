###############################################################################
#
#   Basic manipulation
#
###############################################################################

elem_type{T <: Nemo.RingElem}(::Type{CoefficientRing{T}}) = n_unknown{T}

parent_type{T <: Nemo.RingElem}(::Type{n_unknown{T}}) = CoefficientRing{T}

parent(a::n_unknown) = a.parent

base_ring(R::CoefficientRing) = R.base_ring

base_ring(a::n_unknown) = base_ring(parent(a))

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::CoefficientRing{T}){T <: Nemo.RingElem}(a::T)
   return n_unknown(libSingular.number(a), R)
end

function (R::CoefficientRing{T}){T <: Nemo.RingElem}(a::Integer)
   return R(base_ring(R)(a))
end

function (R::CoefficientRing{T}){T <: Nemo.RingElem}(a::n_unknown{T})
   base_ring(a) != base_ring(R) && error("Unable to coerce element")
   return a
end

function (R::CoefficientRing{T}){T <: Nemo.RingElem, S <: Nemo.RingElem}(a::S)
   U = base_ring(R)
   parent(a) != U && error("Unable to coerce element")
   return R(U(a))
end

###############################################################################
#
#   CoefficientRing constructor
#
###############################################################################

function CoefficientRing(R::Nemo.Ring)
   T = elem_type(R)
   return CoefficientRing{T}(R)
end


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

function check_parent(a::n_unknown, b::n_unknown)
   parent(a) != parent(b) && error("Incompatible parents")
end

function canonical_unit(a::n_unknown)
   R = parent(a)
   n = libSingular.julia(a.ptr)
   return R(libSingular.number(canonical_unit(n)))
end

function deepcopy_internal(a::n_unknown, dict::IdDict)
   R = parent(a)
   n = libSingular.julia(a.ptr)
   return R(libSingular.number(deepcopy(n)))
end

function one(R::CoefficientRing)
   return R(libSingular.n_Init(1, R.ptr))
end

function zero(R::CoefficientRing)
   return R(libSingular.n_Init(0, R.ptr))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::CoefficientRing)
   print(io, "Singular coefficient ring over ")
   show(io, R.base_ring)
end

function show(io::IO, a::n_unknown)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(a.ptr)	
   libSingular.n_Write(nn, parent(a).ptr, false)

   m = libSingular.StringEndS()

   print(io, m)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::n_unknown)
   R = parent(a)
   n = libSingular.n_Neg(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +{T <: Nemo.RingElem}(a::n_unknown{T}, b::n_unknown{T})
   check_parent(a, b)
   R = parent(a)
   n = libSingular.n_Add(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function -{T <: Nemo.RingElem}(a::n_unknown{T}, b::n_unknown{T})
   check_parent(a, b)
   R = parent(a)
   n = libSingular.n_Sub(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function *{T <: Nemo.RingElem}(a::n_unknown{T}, b::n_unknown{T})
   check_parent(a, b)
   R = parent(a)
   n = libSingular.n_Mult(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function =={T <: Nemo.RingElem}(a::n_unknown{T}, b::n_unknown{T})
   check_parent(a, b)
   R = parent(a)
   return libSingular.n_Equal(a.ptr, b.ptr, R.ptr)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv{T <: Nemo.FieldElem}(a::n_unknown{T})
   R = parent(a)
   n = libSingular.n_Invers(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact{T <: Nemo.RingElem}(a::n_unknown{T}, b::n_unknown{T})
   check_parent(a, b)
   R = parent(a)
   n = libSingular.n_ExactDiv(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd{T <: Nemo.RingElem}(a::n_unknown{T}, b::n_unknown{T})
   check_parent(a, b)
   R = parent(a)
   n = libSingular.n_Gcd(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx{T <: Nemo.RingElem}(a::n_unknown{T}, b::n_unknown{T})
   check_parent(a, b)
   R = parent(a)
   s = [libSingular.n_Init(0, R.ptr)]
   t = [libSingular.n_Init(0, R.ptr)]
   n = libSingular.n_ExtGcd(a.ptr, b.ptr, pointer(s), pointer(t), R.ptr)
   return R(n), R(s[]), R(t[])
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule{S <: Nemo.RingElem, T <: Integer}(C::Type{n_unknown{S}}, ::Type{T}) = n_unknown{S}

promote_rule{T <: Nemo.RingElem}(::Type{n_unknown{T}}, ::Type{T}) = n_unknown{T}

promote_rule1{T <: Nemo.RingElem, U <: Nemo.RingElem}(::Type{U}, ::Type{n_unknown{T}}) = promote_rule(U, n_unknown{T})

function promote_rule1{T <: Nemo.RingElem, U <: Nemo.RingElem}(::Type{n_unknown{T}}, ::Type{n_unknown{U}})
   promote_rule(T, n_unknown{U}) == T ? n_unknown{T} : Union{}
end

function promote_rule{T <: Nemo.RingElem, U <: Nemo.RingElem}(::Type{n_unknown{T}}, ::Type{U})
   promote_rule(T, U) == T ? n_unknown{T} : promote_rule1(U, n_unknown{T})
end

promote_rule{S <: Nemo.RingElem}(::Type{n_unknown{S}}, ::Type{Nemo.fmpz}) = n_unknown{S}

promote_rule{S <: Nemo.FieldElem}(::Type{n_unknown{S}}, ::Type{Nemo.fmpq}) = n_unknown{S}

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

function (R::CoefficientRing{T}){T <: Nemo.RingElem}(ptr::libSingular.number)
   return n_unknown{T}(ptr, R)
end

function (R::CoefficientRing{T}){T <: Nemo.RingElem, S <: Nemo.RingElem}(a::S)
   U = base_ring(R)
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


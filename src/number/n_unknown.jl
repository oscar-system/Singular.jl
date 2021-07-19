###############################################################################
#
#   Basic manipulation
#
###############################################################################

elem_type(::Type{CoefficientRing{T}}) where T <: Nemo.RingElem = n_unknown{T}

parent_type(::Type{n_unknown{T}}) where T <: Nemo.RingElem = CoefficientRing{T}

parent(a::n_unknown) = a.parent

base_ring(R::CoefficientRing) = R.base_ring

base_ring(a::n_unknown) = base_ring(parent(a))

function check_parent(a::n_unknown, b::n_unknown)
   parent(a) != parent(b) && error("Incompatible parents")
end

function canonical_unit(a::n_unknown)
   R = parent(a)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return R(canonical_unit(n))
end

function deepcopy_internal(a::n_unknown, dict::IdDict)
   R = parent(a)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return R(deepcopy(n))
end

function one(R::CoefficientRing)
   return R(libSingular.n_Init(1, R.ptr))
end

function zero(R::CoefficientRing)
   return R(libSingular.n_Init(0, R.ptr))
end

function hash(a::n_unknown, h::UInt)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return xor(hash(n, h), 0x664e59de562461fe%UInt)
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
   R = parent(a)
   GC.@preserve a R libSingular.n_Write(a.ptr, R.ptr, false)
   print(io, libSingular.StringEndS())
end

function AbstractAlgebra.expressify(a::n_unknown; context = nothing)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return AbstractAlgebra.expressify(n; context = context)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::n_unknown)
   R = parent(a)
   n = GC.@preserve a R libSingular.n_Neg(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Add(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function -(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Sub(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function *(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Mult(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R return libSingular.n_Equal(a.ptr, b.ptr, R.ptr)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::n_unknown{T}) where T <: Nemo.FieldElem
   R = parent(a)
   n = GC.@preserve a R libSingular.n_Invers(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_ExactDiv(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Gcd(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx(x::n_unknown{T}, y::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R begin
      s = Ref(libSingular.n_Init(0, R.ptr))
      t = Ref(libSingular.n_Init(0, R.ptr))
      g = libSingular.n_ExtGcd(x.ptr, y.ptr, s, t, R.ptr)
      return R(g), R(s[]), R(t[])
   end
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule(C::Type{n_unknown{S}}, ::Type{T}) where {S <: Nemo.RingElem, T <: Integer} = n_unknown{S}

promote_rule(::Type{n_unknown{T}}, ::Type{T}) where {T <: Nemo.RingElem} = n_unknown{T}

promote_rule1(::Type{U}, ::Type{n_unknown{T}}) where {T <: Nemo.RingElem, U <: Nemo.RingElem} = promote_rule(U, n_unknown{T})

function promote_rule1(::Type{n_unknown{T}}, ::Type{n_unknown{U}}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, n_unknown{U}) == T ? n_unknown{T} : Union{}
end

function promote_rule(::Type{n_unknown{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? n_unknown{T} : promote_rule1(U, n_unknown{T})
end

promote_rule(::Type{n_unknown{S}}, ::Type{Nemo.fmpz}) where {S <: Nemo.RingElem} = n_unknown{S}

promote_rule(::Type{n_unknown{S}}, ::Type{Nemo.fmpq}) where {S <: Nemo.FieldElem} = n_unknown{S}

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::CoefficientRing{T})(a::T) where T <: Nemo.RingElem
   return n_unknown(libSingular.cast_void_to_number(libSingular.number(a)), R)
end

function (R::CoefficientRing{T})(a::Integer) where T <: Nemo.RingElem
   return R(base_ring(R)(a))
end

function (R::CoefficientRing{T})(a::Rational) where T <: Nemo.RingElem
   return R(base_ring(R)(a))
end

function (R::CoefficientRing{T})(a::n_Z) where T <: Nemo.RingElem
   return R(base_ring(R)(convert(BigInt, a)))
end

function (R::CoefficientRing{T})(a::n_Q) where T <: Nemo.RingElem
   return R(base_ring(R)(Rational{BigInt}(a)))
end

function (R::CoefficientRing{T})(a::n_unknown{T}) where T <: Nemo.RingElem
   base_ring(a) != base_ring(R) && error("Unable to coerce element")
   return a
end

function (R::CoefficientRing{T})(ptr::libSingular.number_ptr) where {T <: Nemo.RingElem}
   return n_unknown{T}(ptr, R)
end

function (R::CoefficientRing{T})(a::S) where {T <: Nemo.RingElem, S <: Nemo.RingElem}
   U = base_ring(R)
   return R(U(a))
end

###############################################################################
#
#   CoefficientRing constructor
#
###############################################################################

mutable struct MutableRingWrapper{S, T} <: Nemo.Ring
   data::S
end

mutable struct MutableRingElemWrapper{S, T} <: Nemo.RingElem
   data::T
   parent::MutableRingWrapper{S, T}
end

function promote_rule(::Type{n_unknown{MutableRingElemWrapper{S, T}}}, ::Type{T}) where {T <: Nemo.RingElem, S}
   return n_unknown{MutableRingElemWrapper{S, T}}
end

function expressify(a::MutableRingWrapper{S, T}; context = nothing) where {S, T}
   return expressify(a.data, context = context)
end

function expressify(a::MutableRingElemWrapper{S, T}; context = nothing) where {S, T}
   return expressify(a.data, context = context)
end

function Base.show(io::IO, a::MutableRingWrapper{S, T}) where {S, T}
   Base.show(io, a.data)
end

function Base.show(io::IO, a::MutableRingElemWrapper{S, T}) where {S, T}
   Base.show(io, a.data)
end

function elem_type(::Type{MutableRingWrapper{S, T}}) where {S, T}
   return MutableRingElemWrapper{S, T}
end

function parent(a::MutableRingElemWrapper{S, T}) where {S, T}
   return a.parent
end

function (R::MutableRingWrapper{S, T})() where {S, T}
   return MutableRingElemWrapper{S, T}(R.data(), R)
end

function (R::MutableRingWrapper{S, T})(a::MutableRingElemWrapper) where {S, T}
   R == a.parent || error("Unable to coerce element")
   return MutableRingElemWrapper{S, T}(R.data(a.data), R)
end

function Base.deepcopy_internal(a::MutableRingElemWrapper{S, T}, dict::IdDict) where {S, T}
   return MutableRingElemWrapper{S, T}(deepcopy_internal(a.data, dict), a.parent)
end

# should be R::S
function (R::Nemo.Ring)(a::n_unknown{MutableRingElemWrapper{S, T}}) where {S, T}
   GC.@preserve a begin
      ja = libSingular.julia(libSingular.cast_number_to_void(a.ptr))
      R == ja.parent.data || error("Unable to coerce element")
      return ja.data::T
   end
end

function (R::MutableRingWrapper{S, T})(a) where {S, T}
   return MutableRingElemWrapper{S, T}(R.data(a), R)
end

function zero(R::MutableRingWrapper{S, T}) where {S, T}
   return MutableRingElemWrapper{S, T}(zero(R.data), R)
end

function one(R::MutableRingWrapper{S, T}) where {S, T}
   return MutableRingElemWrapper{S, T}(one(R.data), R)
end

function Base.hash(a::MutableRingElemWrapper{S, T}, b::UInt) where {S, T}
   return hash(a.data, b)
end

# one input, one non-wrapped output
for op in (:iszero, :isone)
   @eval begin
      function ($op)(a::MutableRingElemWrapper{S, T}) where {S, T}
         return ($op)(a.data)
      end
   end
end

# one input, one wrapped output
for op in (:-, :zero, :one)
   @eval begin
      function ($op)(a::MutableRingElemWrapper{S, T}) where {S, T}
         return MutableRingElemWrapper{S, T}(($op)(a.data), a.parent)
      end
   end
end

# two inputs, one non-wrapped output
for op in (:(==), )
   @eval begin
      function ($op)(a::MutableRingElemWrapper{S, T}, b::MutableRingElemWrapper{S, T}) where {S, T}
         return ($op)(a.data, b.data)
      end
   end
end

# two inputs, one wrapped output
for op in (:+, :-, :*, :div, :divexact, :gcd)
   @eval begin
      function ($op)(a::MutableRingElemWrapper{S, T}, b::MutableRingElemWrapper{S, T}) where {S, T}
         return MutableRingElemWrapper{S, T}(($op)(a.data, b.data), a.parent)
      end
   end
end

# two inputs, one non-wrapped and one wrapped output
for op in (:divides, )
   @eval begin
      function ($op)(a::MutableRingElemWrapper{S, T}, b::MutableRingElemWrapper{S, T}) where {S, T}
         res1, res2 = ($op)(a.data, b.data)
         return (res1, MutableRingElemWrapper{S, T}(res2, a.parent))
      end
   end
end

# two inputs, two wrapped outputs
for op in (:divrem, )
   @eval begin
      function ($op)(a::MutableRingElemWrapper{S, T}, b::MutableRingElemWrapper{S, T}) where {S, T}
         res1, res2 = ($op)(a.data, b.data)
         return (MutableRingElemWrapper{S, T}(res1, a.parent),
                 MutableRingElemWrapper{S, T}(res2, a.parent))
      end
   end
end

# two inputs, three wrapped outputs
for op in (:gcdx, )
   @eval begin
      function ($op)(a::MutableRingElemWrapper{S, T}, b::MutableRingElemWrapper{S, T}) where {S, T}
         res1, res2, res3 = ($op)(a.data, b.data)
         return (MutableRingElemWrapper{S, T}(res1, a.parent),
                 MutableRingElemWrapper{S, T}(res2, a.parent),
                 MutableRingElemWrapper{S, T}(res3, a.parent))
      end
   end
end

function addeq!(z::MutableRingElemWrapper{S, T}, a::MutableRingElemWrapper{S, T}) where {S, T}
   addeq!(z.data, a.data)
   return z
end

function mul!(z::MutableRingElemWrapper{S, T}, a::MutableRingElemWrapper{S, T}, b::MutableRingElemWrapper{S, T}) where {S, T}
   mul!(z.data, a.data, b.data)
   return z
end


function CoefficientRing(R::Nemo.Ring)
   T = elem_type(R)

   if VERSION >= v"1.8"
      ok = ismutabletype(T)
   elseif VERSION >= v"1.5"
      ok = ismutable(R())
   else
      ok = !isimmutable(R())
   end

   if ok
      return CoefficientRing{T}(R)
   else
      RR = MutableRingWrapper{typeof(R), elem_type(R)}(R)
      return CoefficientRing{elem_type(RR)}(RR)
   end
end


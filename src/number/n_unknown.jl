###############################################################################
#
#   Basic manipulation
#
###############################################################################
elem_type(::Type{N_Ring{T}}) where T <: Nemo.RingElem = n_RingElem{T}
elem_type(::Type{N_Field{T}}) where T <: Nemo.FieldElem = n_FieldElem{T}

parent_type(::Type{n_RingElem{T}}) where T <: Nemo.RingElem = N_Ring{T}
parent_type(::Type{n_FieldElem{T}}) where T <: Nemo.RingElem = N_Field{T}

parent(a::n_unknown) = a.parent

base_ring(R::N_unknown) = R.base_ring

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

function one(R::N_unknown)
   return R(libSingular.n_Init(1, R.ptr))
end

function zero(R::N_unknown)
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

function show(io::IO, R::N_unknown)
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

function +(a::n_unknown, b::n_unknown)
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Add(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function -(a::n_unknown, b::n_unknown)
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Sub(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function *(a::n_unknown, b::n_unknown)
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

function ==(a::n_unknown, b::n_unknown)
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R return libSingular.n_Equal(a.ptr, b.ptr, R.ptr)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::n_unknown)
   R = parent(a)
   n = GC.@preserve a R libSingular.n_Invers(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::n_unknown, b::n_unknown)
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

function gcd(a::n_unknown, b::n_unknown)
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

function gcdx(x::n_unknown, y::n_unknown)
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

promote_rule(C::Type{n_RingElem{S}}, ::Type{T}) where {S <: Nemo.RingElem, T <: Integer} = n_RingElem{S}
promote_rule(C::Type{n_FieldElem{S}}, ::Type{T}) where {S <: Nemo.FieldElem, T <: Integer} = n_FieldElem{S}

promote_rule(::Type{n_RingElem{T}}, ::Type{T}) where {T <: Nemo.RingElem} = n_RingElem{T}
promote_rule(::Type{n_FieldElem{T}}, ::Type{T}) where {T <: Nemo.FieldElem} = n_FieldElem{T}

promote_rule1(::Type{U}, ::Type{n_RingElem{T}}) where {T <: Nemo.RingElem, U <: Nemo.RingElem} = promote_rule(U, n_RingElem{T})
promote_rule1(::Type{U}, ::Type{n_FieldElem{T}}) where {T <: Nemo.FieldElem, U <: Nemo.FieldElem} = promote_rule(U, n_FieldElem{T})

function promote_rule1(::Type{n_RingElem{T}}, ::Type{n_RingElem{U}}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, n_RingElem{U}) == T ? n_RingElem{T} : Union{}
end

function promote_rule1(::Type{n_FieldElem{T}}, ::Type{n_FieldElem{U}}) where {T <: Nemo.FieldElem, U <: Nemo.FieldElem}
   promote_rule(T, n_FieldElem{U}) == T ? n_FieldElem{T} : Union{}
end

function promote_rule(::Type{n_RingElem{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? n_RingElem{T} : promote_rule1(U, n_RingElem{T})
end

function promote_rule(::Type{n_FieldElem{T}}, ::Type{U}) where {T <: Nemo.FieldElem, U <: Nemo.FieldElem}
   promote_rule(T, U) == T ? n_FieldElem{T} : promote_rule1(U, n_FieldElem{T})
end

promote_rule(::Type{n_RingElem{S}}, ::Type{Nemo.fmpz}) where {S <: Nemo.RingElem} = n_RingElem{S}
promote_rule(::Type{n_FieldElem{S}}, ::Type{Nemo.fmpz}) where {S <: Nemo.FieldElem} = n_FieldElem{S}

promote_rule(::Type{n_RingElem{S}}, ::Type{Nemo.fmpq}) where {S <: Nemo.RingElem} = n_RingElem{S}
promote_rule(::Type{n_FieldElem{S}}, ::Type{Nemo.fmpq}) where {S <: Nemo.FieldElem} = n_FieldElem{S}

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::N_Ring{T})(a::T) where T <: Nemo.RingElem
   return n_RingElem(libSingular.cast_void_to_number(libSingular.number(a)), R)
end

function (R::N_Field{T})(a::T) where T <: Nemo.FieldElem
   return n_FieldElem(libSingular.cast_void_to_number(libSingular.number(a)), R)
end

function (R::N_Ring{T})(a::Integer) where T <: Nemo.RingElem
   return R(base_ring(R)(a))
end

function (R::N_Field{T})(a::Integer) where T <: Nemo.FieldElem
   return R(base_ring(R)(a))
end

function (R::N_Ring{T})(a::Rational) where T <: Nemo.RingElem
   return R(base_ring(R)(a))
end

function (R::N_Field{T})(a::Rational) where T <: Nemo.FieldElem
   return R(base_ring(R)(a))
end

function (R::N_Ring{T})(a::n_Z) where T <: Nemo.RingElem
   return R(base_ring(R)(convert(BigInt, a)))
end

function (R::N_Field{T})(a::n_Z) where T <: Nemo.FieldElem
   return R(base_ring(R)(convert(BigInt, a)))
end

function (R::N_Ring{T})(a::n_Q) where T <: Nemo.RingElem
   return R(base_ring(R)(Rational{BigInt}(a)))
end

function (R::N_Field{T})(a::n_Q) where T <: Nemo.FieldElem
   return R(base_ring(R)(Rational{BigInt}(a)))
end

function (R::N_Ring{T})(a::n_RingElem{T}) where T <: Nemo.RingElem
   base_ring(a) != base_ring(R) && error("Unable to coerce element")
   return a
end

function (R::N_Field{T})(a::n_FieldElem{T}) where T <: Nemo.RingElem
   base_ring(a) != base_ring(R) && error("Unable to coerce element")
   return a
end

function (R::N_Ring{T})(ptr::libSingular.number_ptr) where T <: Nemo.RingElem
   return n_RingElem{T}(ptr, R)
end

function (R::N_Field{T})(ptr::libSingular.number_ptr) where T <: Nemo.FieldElem
   return n_FieldElem{T}(ptr, R)
end

function (R::N_Ring{T})(a::S) where {T <: Nemo.RingElem, S <: Nemo.RingElem}
   U = base_ring(R)
   return R(U(a))
end

function (R::N_Field{T})(a::S) where {T <: Nemo.FieldElem, S <: Nemo.RingElem}
   U = base_ring(R)
   return R(U(a))
end

###############################################################################
#
#   N_Ring N_Field constructor
#
###############################################################################

# mutable
mutable struct RingWrapper{S, T} <: Nemo.Ring
   data::S
end

mutable struct FieldWrapper{S, T} <: Nemo.Field
   data::S
end

mutable struct RingElemWrapper{S, T} <: Nemo.RingElem
   data::T
   parent::RingWrapper{S, T}
end

mutable struct FieldElemWrapper{S, T} <: Nemo.FieldElem
   data::T
   parent::FieldWrapper{S, T}
end

function promote_rule(::Type{n_RingElem{RingElemWrapper{S, T}}}, ::Type{T}) where {T <: Nemo.RingElem, S}
   return n_RingElem{RingElemWrapper{S, T}}
end

function promote_rule(::Type{n_FieldElem{FieldElemWrapper{S, T}}}, ::Type{T}) where {T <: Nemo.FieldElem, S}
   return n_FieldElem{FieldElemWrapper{S, T}}
end

function expressify(a::Union{RingWrapper, FieldWrapper}; context = nothing)
   return expressify(a.data, context = context)
end

function expressify(a::Union{RingElemWrapper, FieldElemWrapper}; context = nothing)
   return expressify(a.data, context = context)
end

function Base.show(io::IO, a::Union{RingWrapper, FieldWrapper})
   Base.show(io, a.data)
end

function Base.show(io::IO, a::Union{RingElemWrapper, FieldElemWrapper})
   Base.show(io, a.data)
end

elem_type(::Type{RingWrapper{S, T}}) where {S, T} = RingElemWrapper{S, T}
elem_type(::Type{FieldWrapper{S, T}}) where {S, T} = FieldElemWrapper{S, T}

parent(a::RingElemWrapper{S, T}) where {S, T} = a.parent
parent(a::FieldElemWrapper{S, T}) where {S, T} = a.parent

function (R::RingWrapper{S, T})() where {S, T}
   return RingElemWrapper{S, T}(R.data(), R)
end

function (R::FieldWrapper{S, T})() where {S, T}
   return FieldElemWrapper{S, T}(R.data(), R)
end

function (R::RingWrapper{S, T})(a::RingElemWrapper) where {S, T}
   R == a.parent || error("Unable to coerce element")
   return RingElemWrapper{S, T}(R.data(a.data), R)
end

function (R::FieldWrapper{S, T})(a::FieldElemWrapper) where {S, T}
   R == a.parent || error("Unable to coerce element")
   return FieldElemWrapper{S, T}(R.data(a.data), R)
end

# should be R::S
function (R::Nemo.Ring)(a::n_RingElem{RingElemWrapper{S, T}}) where {S, T}
   GC.@preserve a begin
      ja = libSingular.julia(libSingular.cast_number_to_void(a.ptr))
      R == ja.parent.data || error("Unable to coerce element")
      return ja.data::T
   end
end

function (R::Nemo.Field)(a::n_FieldElem{FieldElemWrapper{S, T}}) where {S, T}
   GC.@preserve a begin
      ja = libSingular.julia(libSingular.cast_number_to_void(a.ptr))
      R == ja.parent.data || error("Unable to coerce element")
      return ja.data::T
   end
end

for (rw, rew) in [(:RingWrapper, :RingElemWrapper),
                  (:FieldWrapper, :FieldElemWrapper)]

@eval begin
   function Base.deepcopy_internal(a::($rew){S, T}, dict::IdDict) where {S, T}
      return ($rew){S, T}(deepcopy_internal(a.data, dict), a.parent)
   end

   function Base.hash(a::($rew){S, T}, b::UInt) where {S, T}
      return hash(a.data, b)
   end

   function (R::($rw){S, T})(a) where {S, T}
      return ($rew){S, T}(R.data(a), R)
   end

   function zero(R::($rw){S, T}) where {S, T}
      return ($rew){S, T}(zero(R.data), R)
   end

   function one(R::($rw){S, T}) where {S, T}
      return ($rew){S, T}(one(R.data), R)
   end

   function addeq!(z::($rew){S, T}, a::($rew){S, T}) where {S, T}
      addeq!(z.data, a.data)
      return z
   end

   function mul!(z::($rew){S, T}, a::($rew){S, T}, b::($rew){S, T}) where {S, T}
      mul!(z.data, a.data, b.data)
      return z
   end
end

# one input, one non-wrapped output
for op in (:iszero, :isone)
   @eval begin
      function ($op)(a::($rew){S, T}) where {S, T}
         return ($op)(a.data)
      end
   end
end

# one input, one wrapped output
for op in (:-, :zero, :one)
   @eval begin
      function ($op)(a::($rew){S, T}) where {S, T}
         return ($rew){S, T}(($op)(a.data), a.parent)
      end
   end
end

# two inputs, one non-wrapped output
for op in (:(==), )
   @eval begin
      function ($op)(a::($rew){S, T}, b::($rew){S, T}) where {S, T}
         return ($op)(a.data, b.data)
      end
   end
end

# two inputs, one wrapped output
for op in (:+, :-, :*, :div, :divexact, :gcd)
   @eval begin
      function ($op)(a::($rew){S, T}, b::($rew){S, T}) where {S, T}
         return ($rew){S, T}(($op)(a.data, b.data), a.parent)
      end
   end
end

# two inputs, one non-wrapped and one wrapped output
for op in (:divides, )
   @eval begin
      function ($op)(a::($rew){S, T}, b::($rew){S, T}) where {S, T}
         res1, res2 = ($op)(a.data, b.data)
         return (res1, ($rew){S, T}(res2, a.parent))
      end
   end
end

# two inputs, two wrapped outputs
for op in (:divrem, )
   @eval begin
      function ($op)(a::($rew){S, T}, b::($rew){S, T}) where {S, T}
         res1, res2 = ($op)(a.data, b.data)
         return (($rew){S, T}(res1, a.parent),
                 ($rew){S, T}(res2, a.parent))
      end
   end
end

# two inputs, three wrapped outputs
for op in (:gcdx, )
   @eval begin
      function ($op)(a::($rew){S, T}, b::($rew){S, T}) where {S, T}
         res1, res2, res3 = ($op)(a.data, b.data)
         return (($rew){S, T}(res1, a.parent),
                 ($rew){S, T}(res2, a.parent),
                 ($rew){S, T}(res3, a.parent))
      end
   end
end

end # for (rw, rew) in


function CoefficientRing(R::Nemo.Ring)
   T = elem_type(R)

   if VERSION >= v"1.8"
      ok = ismutabletype(T)
   elseif VERSION >= v"1.5"
      ok = ismutable(R())
   else
      ok = !isimmutable(R())
   end

   if R isa Nemo.Field
      if ok
         return N_Field{T}(R)
      else
         RR = FieldWrapper{typeof(R), elem_type(R)}(R)
         return N_Field{elem_type(RR)}(RR)
      end
   else
      if ok
         return N_Ring{T}(R)
      else
         RR = RingWrapper{typeof(R), elem_type(R)}(R)
         return N_Ring{elem_type(RR)}(RR)
      end
   end
end

export ngens, coeff, isgen, content, primpart, lead_exponent

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(p::spoly) = p.parent

base_ring(R::PolyRing) = R.base_ring

base_ring(p::spoly) = base_ring(parent(p))

elem_type{T <: Nemo.RingElem}(::PolyRing{T}) = spoly{T}

parent_type{T <: Nemo.RingElem}(::Type{spoly{T}}) = PolyRing{T}

ngens(R::PolyRing) = Int(libSingular.rVar(R.ptr))

characteristic(R::PolyRing) = Int(libSingular.rChar(R.ptr))

function gens(R::PolyRing)
   n = ngens(R)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

zero(R::PolyRing) = R()

one(R::PolyRing) = R(1)

iszero(p::spoly) = p.ptr == C_NULL

function isone(p::spoly)
   return Bool(libSingular.p_IsOne(p.ptr, parent(p).ptr))
end

function isgen(p::spoly)
   R = parent(p)
   if p.ptr == C_NULL || libSingular.pNext(p.ptr) != C_NULL || 
     !Bool(libSingular.n_IsOne(libSingular.pGetCoeff(p.ptr), base_ring(p).ptr))
      return false
   end
   n = 0
   for i = 1:ngens(R)
      d = libSingular.p_GetExp(p.ptr, Cint(i), R.ptr)
      if d > 1
         return false
      elseif d == 1
         if n == 1
            return false
         end
         n = 1
      end
   end
   return true   
end

function isconstant(p::spoly)
   R = parent(p)
   if p.ptr == C_NULL
      return true
   end
   if libSingular.pNext(p.ptr) != C_NULL
      return false
   end
   for i = 1:ngens(R)
      if libSingular.p_GetExp(p.ptr, Cint(i), R.ptr) != 0
         return false
      end
   end
   return true   
end

function isunit(p::spoly)
   return p.ptr != C_NULL && libSingular.pNext(p.ptr) == C_NULL &&
          Bool(libSingular.p_IsUnit(p.ptr, parent(p).ptr))
end

length(p::spoly) = Int(libSingular.pLength(p.ptr))

function coeff(p::spoly, i::Int)
   i = length(p) - i - 1
   R = base_ring(p)
   if i < 0 || iszero(p)
      return R()
   end
   ptr = p.ptr
   for i = 1:i
      ptr = libSingular.pNext(ptr)
      if ptr == C_NULL
         return R()
      end
   end
   return R(libSingular.n_Copy(libSingular.pGetCoeff(ptr), R.ptr))
end

function exponent(p::spoly, i::Int)
   i = length(p) - i - 1
   R = parent(p)
   n = ngens(R)
   if i < 0 || iszero(p)
      return zeros(Int, n)
   end
   ptr = p.ptr
   for i = 1:i
      ptr = libSingular.pNext(ptr)
      if ptr == C_NULL
         return zeros(Int, n)
      end
   end
   A = Array{Int}(n)
   libSingular.p_GetExpVL(ptr, A, R.ptr)
   return A
end

function exponent!(A::Array{Int, 1}, p::spoly, i::Int)
   i = length(p) - i - 1
   R = parent(p)
   n = ngens(R)
   @assert length(A) == n
   if i < 0 || iszero(p)
      for i=1:n
        A[i] = 0
      end
      return A
   end
   ptr = p.ptr
   for i = 1:i
      ptr = libSingular.pNext(ptr)
      if ptr == C_NULL
        for i=1:n
          A[i] = 0
        end
      end
   end
   libSingular.p_GetExpVL(ptr, A, R.ptr)
   return A
end

mutable struct coeffs_expos
  E::Array{Int, 1}
  c::Nemo.RingElem
  R::Nemo.Ring
  Rx::Nemo.Ring
  a::spoly
  function coeffs_expos(p::spoly)
    r = new()
    Rx = parent(p)
    n = ngens(Rx)
    r.E = zeros(Int, n)
    r.Rx = Rx
    r.R = base_ring(p)
    r.c = r.R(0)
    r.a = p
    return r
  end
end

function show(io::IO, sp::coeffs_expos)
  println(io, "Coefficients and exponent iterator for $(sp.a)")
end

function Base.start(sp::coeffs_expos)
  return sp.a.ptr
end

function Base.next(sp::coeffs_expos, p)
  if p == C_NULL
    return (sp.c, sp.E), p
  end

  libSingular.p_GetExpVL(p, sp.E, sp.Rx.ptr)
  sp.c = sp.R(libSingular.n_Copy(libSingular.pGetCoeff(p), sp.R.ptr))

  return (sp.c, sp.E), libSingular.pNext(p) 
end

function Base.done(sp::coeffs_expos, p)
  return p == C_NULL
end

function lead_exponent(p::spoly)
   R = parent(p)
   n = ngens(R)
   A = Array{Int}(n)
   libSingular.p_GetExpVL(p.ptr, A, R.ptr)
   return A
end

function deepcopy(p::spoly)
   p2 = libSingular.p_Copy(p.ptr, parent(p).ptr)
   return parent(p)(p2)
end

function check_parent{T <: Nemo.RingElem}(a::spoly{T}, b::spoly{T})
   parent(a) != parent(b) && error("Incompatible parent objects")
end

function canonical_unit{T <: Nemo.RingElem}(a::spoly{T})
  return a == 0 ? one(base_ring(a)) : canonical_unit(coeff(a, 0))
end
   
###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::PolyRing)
   m = libSingular.rString(R.ptr)
   s = unsafe_string(m)
   libSingular.omFree(Ptr{Void}(m))
   print(io, "Singular Polynomial Ring ", s)
end

function show(io::IO, a::spoly)
   m = libSingular.p_String(a.ptr, parent(a).ptr)
   s = unsafe_string(m)
   libSingular.omFree(Ptr{Void}(m))
   print(io, s)
end

show_minus_one{T <: Nemo.RingElem}(::Type{spoly{T}}) = show_minus_one(T)

needs_parentheses(x::spoly) = length(x) > 1

isnegative(x::spoly) = isconstant(x) && !iszero(x) && isnegative(coeff(x, 0))

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(a::spoly)
   a1 = libSingular.p_Copy(a.ptr, parent(a).ptr)
   s = libSingular.p_Neg(a1, parent(a).ptr)
   return parent(a)(s) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +{T <: Nemo.RingElem}(a::spoly{T}, b::spoly{T})
   check_parent(a, b)
   a1 = libSingular.p_Copy(a.ptr, parent(a).ptr)
   b1 = libSingular.p_Copy(b.ptr, parent(a).ptr)
   s = libSingular.p_Add_q(a1, b1, parent(a).ptr)
   return parent(a)(s) 
end

function -{T <: Nemo.RingElem}(a::spoly{T}, b::spoly{T})
   check_parent(a, b)
   a1 = libSingular.p_Copy(a.ptr, parent(a).ptr)
   b1 = libSingular.p_Copy(b.ptr, parent(a).ptr)
   s = libSingular.p_Sub(a1, b1, parent(a).ptr)
   return parent(a)(s) 
end

function *{T <: Nemo.RingElem}(a::spoly{T}, b::spoly{T})
   check_parent(a, b)
   a1 = libSingular.p_Copy(a.ptr, parent(a).ptr)
   b1 = libSingular.p_Copy(b.ptr, parent(a).ptr)
   s = libSingular.p_Mult_q(a1, b1, parent(a).ptr)
   return parent(a)(s) 
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::spoly, y::Int)
    y > typemax(Cint) || y < typemin(Cint) && throw(DomainError())
    R = parent(x)
    if isone(x)
       return deepcopy(x)
    elseif y == 0
       return one(R)
    elseif y == 1
       return deepcopy(x)
    end
    x1 = libSingular.p_Copy(x.ptr, R.ptr)
    p = libSingular.p_Power(x1, Cint(y), R.ptr)
    return R(p)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function =={T <: Nemo.RingElem}(x::spoly{T}, y::spoly{T})
    check_parent(x, y)
    return Bool(libSingular.p_EqualPolys(x.ptr, y.ptr, parent(x).ptr))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::spoly, y::spoly)
   check_parent(x, y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   p = libSingular.p_Divide(x1, y1, R.ptr)
   return R(p)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd{T <: Nemo.RingElem}(x::spoly{T}, y::spoly{T})
   check_parent(x, y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   p = libSingular.singclap_gcd(x1, y1, R.ptr)
   return R(p)
end

function gcdx{T <: Nemo.FieldElem}(x::spoly{T}, y::spoly{T})
   check_parent(x, y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   s = libSingular.poly_ref(libSingular.poly(C_NULL))
   t = libSingular.poly_ref(libSingular.poly(C_NULL))
   p = libSingular.poly_ref(libSingular.poly(C_NULL))
   libSingular.singclap_extgcd(x1, y1, p, s, t, R.ptr)
   return R(p[]), R(s[]), R(t[])
end

function lcm{T <: Nemo.RingElem}(x::spoly{T}, y::spoly{T})
   if iszero(x) && iszero(y)
      return parent(x)()
   end
   return divexact(x*y, gcd(x, y))
end

function primpart(x::spoly)
   R = parent(x)
   p = deepcopy(x)
   libSingular.p_Content(p.ptr, R.ptr)
   return p
end

function content(x::spoly)
   R = base_ring(x)
   d = R()
   for i = 1:length(x)
      d = gcd(d, coeff(x, i - 1))
      if isone(d)
         break
      end
   end
   return d
end

###############################################################################
#
#   Unsafe operations
#
###############################################################################

function addeq!(x::spoly, y::spoly)
    R = parent(x)
    if y.ptr == C_NULL
    elseif x.ptr == C_NULL
       x.ptr = libSingular.p_Copy(y.ptr, R.ptr)
    else
       x.ptr = libSingular.p_Add_q(x.ptr, libSingular.p_Copy(y.ptr, R.ptr), R.ptr)
    end
    return x
end

function mul!(c::spoly, x::spoly, y::spoly) 
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   ptr = libSingular.p_Mult_q(x1, y1, R.ptr)
   if c.ptr != C_NULL
      libSingular.p_Delete(c.ptr, R.ptr)
   end
   c.ptr = ptr
   return c
end

function add!(c::spoly, x::spoly, y::spoly) 
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   ptr = libSingular.p_Add_q(x1, y1, R.ptr)
   if c.ptr != C_NULL
      libSingular.p_Delete(c.ptr, R.ptr)
   end
   c.ptr = ptr
   return c
end

function zero!(x::spoly)
   if x.ptr != C_NULL
      libSingular.p_Delete(x.ptr, parent(x).ptr)
      x.ptr = libSingular.p_ISet(0, parent(x).ptr)
   end
   return x
end

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{spoly{T}}, ::Type{spoly{T}}) where T <: Nemo.RingElem = spoly{T}

function promote_rule(::Type{spoly{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? spoly{T} : Union{}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::PolyRing)()
   T = elem_type(base_ring(R))
   return spoly{T}(R)
end

function (R::PolyRing)(n::Int)
   T = elem_type(base_ring(R))
   return spoly{T}(R, n)
end

function (R::PolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   return spoly{T}(R, BigInt(n))
end

function (R::PolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   return spoly{T}(R, ptr)
end

function (R::PolyRing)(n::libSingular.poly)
   T = elem_type(base_ring(R))
   return spoly{T}(R, n)
end

function (R::PolyRing{T}){T <: Nemo.RingElem}(n::T)
   parent(n) != base_ring(R) && error("Unable to coerce into polynomial ring")
   return spoly{T}(R, n.ptr)
end

function (R::PolyRing{S}){S <: Nemo.RingElem, T <: Nemo.RingElem}(n::T)
   return R(base_ring(R)(n))
end

function (R::PolyRing)(p::spoly)
   parent(p) != R && error("Unable to coerce polynomial")
   return p
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::Union{Ring, Field}, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :degrevlex)
   U = [Symbol(v) for v in s]
   T = elem_type(R)
   parent_obj = PolyRing{T}(R, U, cached, sym2ringorder[ordering])
   return tuple(parent_obj, gens(parent_obj))
end

function PolynomialRing(R::Nemo.Ring, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :degrevlex)
   S = CoefficientRing(R)
   U = [Symbol(v) for v in s]
   T = elem_type(S)
   parent_obj = PolyRing{T}(S, U, cached, sym2ringorder[ordering])
   return tuple(parent_obj, gens(parent_obj))
end


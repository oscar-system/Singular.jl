export SingularPolynomialRing, ngens, coeff, isgen, content, primpart

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(p::spoly) = p.parent

base_ring(R::SingularPolyRing) = R.base_ring

base_ring(p::spoly) = base_ring(parent(p))

elem_type{T <: Nemo.RingElem}(::SingularPolyRing{T}) = spoly{T}

parent_type{T <: Nemo.RingElem}(a::spoly{T}) = SingularPolyRing{T}

ngens(R::SingularPolyRing) = Int(libSingular.rVar(R.ptr))

characteristic(R::SingularPolyRing) = Int(libSingular.rChar(R.ptr))

function gens(R::SingularPolyRing)
   n = ngens(R)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

zero(R::SingularPolyRing) = R()

one(R::SingularPolyRing) = R(1)

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

function deepcopy(p::spoly)
   p2 = libSingular.p_Copy(p.ptr, parent(p).ptr)
   return parent(p)(p2)
end

function check_parent{T <: Nemo.RingElem}(a::spoly{T}, b::spoly{T})
   parent(a) != parent(b) && error("Incompatible parent objects")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::SingularPolyRing)
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

is_negative(x::spoly) = isconstant(x) && !iszero(x) && is_negative(coeff(x, 0))

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
   p = libSingular.singclap_pdivide(x1, y1, R.ptr)
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
    nothing
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
   nothing
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
   nothing
end

function zero!(x::spoly)
   if x.ptr != C_NULL
      libSingular.p_Delete(x.ptr, parent(x).ptr)
      x.ptr = libSingular.p_ISet(0, parent(x).ptr)
   end
   nothing
end

###############################################################################
#
#   Promote rules
#
###############################################################################

Base.promote_rule{U <: Integer, T <: Nemo.RingElem}(::Type{spoly{T}}, ::Type{U}) = spoly{T}

Base.promote_rule{T <: Nemo.RingElem}(::Type{spoly{T}}, ::Type{n_Z}) = spoly{T}

Base.promote_rule{T <: Nemo.RingElem}(::Type{spoly{T}}, ::Type{T}) = spoly{T}

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::SingularPolyRing)()
   T = elem_type(base_ring(R))
   return spoly{T}(R)
end

function (R::SingularPolyRing)(n::Int)
   T = elem_type(base_ring(R))
   return spoly{T}(R, n)
end

function (R::SingularPolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   return spoly{T}(R, BigInt(n))
end

function (R::SingularPolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   return spoly{T}(R, ptr)
end

function (R::SingularPolyRing)(n::libSingular.poly)
   T = elem_type(base_ring(R))
   return spoly{T}(R, n)
end

function (R::SingularPolyRing){T <: Nemo.RingElem}(n::T)
   parent(n) != base_ring(R) && error("Unable to coerce into polynomial ring")
   return spoly{T}(R, n.ptr)
end

function (R::SingularPolyRing)(p::spoly)
   parent(p) != R && error("Unable to coerce polynomial")
   return p
end

###############################################################################
#
#   SingularPolynomialRing constructor
#
###############################################################################

function SingularPolynomialRing(R::Nemo.Ring, s::Array{String, 1}; cached::Bool = true, ordering::Symbol = :degrevlex)
   U = [Symbol(v) for v in s]
   T = elem_type(R)
   parent_obj = SingularPolyRing{T}(R, U, cached, sym2ringorder[ordering])
   return tuple(parent_obj, gens(parent_obj)...)
end


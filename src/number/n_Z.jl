export crt

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::SingularIntegerRing) = n_Z

parent(a::n_Z) = SingularZZ

parent_type(::Type{n_Z}) = SingularIntegerRing

base_ring(a::n_Z) = SingularZZ

base_ring(a::SingularIntegerRing) = SingularZZ

function deepcopy(a::n_Z)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(::SingularIntegerRing) = SingularZZ(1)

zero(::SingularIntegerRing) = SingularZZ(0)

function isone(n::n_Z)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Z)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

isunit(n::n_Z) = !iszero(n)

function num(n::n_Z)
   par = parent(n)
   r = libSingular.n_Copy(n.ptr, par.ptr)
   return par(r)
end

function den(n::n_Z)
   return one(parent(n))
end

function abs(n::n_Z)
   if libSingular.n_GreaterZero(n.ptr, parent(n).ptr) || iszero(n)
      return deepcopy(n)
   else
      return -n
   end
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_Z) = isnegative(x) ? -one(parent(x)) : one(parent(x))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::SingularIntegerRing)
   print(io, "Integer Ring")
end

function show(io::IO, n::n_Z)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
   n.ptr = nn[]

   m = libSingular.StringEndS()
   s = unsafe_string(m) 
   libSingular.omFree(Ptr{Void}(m))

   print(io, s)
end

needs_parentheses(x::n_Z) = false

is_negative(x::n_Z) = !libSingular.n_GreaterZero(x.ptr, parent(x).ptr) &&
                      !iszero(x)

show_minus_one(::Type{n_Z}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Z) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Z, y::n_Z)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Z, y::n_Z)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Z, y::n_Z)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_Z, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_Z) = parent(y)(x) + y

-(x::n_Z, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_Z) = parent(y)(x) - y

*(x::n_Z, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_Z) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::n_Z, y::n_Z)
    libSingular.n_Greater(y.ptr, x.ptr, parent(x).ptr)
end

function ==(x::n_Z, y::n_Z)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end

isequal(x::n_Z, y::n_Z) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_Z, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_Z) = (parent(y)(x) == y)

isequal(x::n_Z, y::Integer) = (x == y)

isequal(x::Integer, y::n_Z) = (x == y)

isless(x::n_Z, y::Integer) = isless(x, parent(x)(y))

isless(x::Integer, y::n_Z) = isless(parent(y)(x), y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_Z, y::Int)
    y < 0 && throw(DomainError())
    if isone(x)
       return x
    elseif y == 0
       return one(parent(x))
    elseif y == 1
       return x
    else
       p = libSingular.n_Power(x.ptr, y, parent(x).ptr)
       return parent(x)(p)
    end
end

###############################################################################
#
#   Exact division
#
###############################################################################

function div(x::n_Z, y::n_Z)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_Z, y::n_Z)
   c = parent(x)
   p = libSingular.n_ExactDiv(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::n_Z, y::n_Z)
   par = parent(x)
   r = [libSingular.n_Init(0, par.ptr)]
   q = libSingular.n_QuotRem(x.ptr, y.ptr, pointer(r), par.ptr)
   return par(q), par(r[])
end

function rem(x::n_Z, y::n_Z)
   par = parent(x)
   r = libSingular.n_IntMod(x.ptr, y.ptr, par.ptr)
   return par(r)
end

mod(x::n_Z, y::n_Z) = rem(x, y)

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_Z, y::n_Z)
   par = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

function lcm(x::n_Z, y::n_Z)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   return div(x*y, gcd(x, y))
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx(x::n_Z, y::n_Z)
   par = parent(x)
   s = [libSingular.n_Init(0, par.ptr)]
   t = [libSingular.n_Init(0, par.ptr)]
   g = libSingular.n_ExtGcd(x.ptr, y.ptr, pointer(s), pointer(t), par.ptr)
   return par(g), par(s[]), par(t[])
end

###############################################################################
#
#   Chinese remainder
#
###############################################################################

function crt(r1::n_Z, m1::n_Z, r2::n_Z, m2::n_Z, signed=false)
   par = parent(r1)
   r = [r1.ptr, r2.ptr]
   m = [m1.ptr, m2.ptr]
   a = libSingular.n_ChineseRemainderSym(r, m, Cint(2), Cint(signed), par.ptr)
   return par(a)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Z, y::n_Z)
    xx = libSingular.number_ref(x.ptr)
    libSingular.n_InpAdd(xx, y.ptr, parent(x).ptr)
    x.ptr = xx[]
    nothing
end

function mul!(x::n_Z, y::n_Z, z::n_Z)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function add!(x::n_Z, y::n_Z, z::n_Z)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function zero!(x::n_Z)
   ptr = libSingular.n_Init(0, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(C::Type{n_Z}, ::Type{T}) = n_Z

###############################################################################
#
#   Parent call functions
#
###############################################################################

(::SingularIntegerRing)() = n_Z()

(R::SingularIntegerRing)(x::Integer) = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 

(::SingularIntegerRing)(n::Int) = n_Z(n)

(::SingularIntegerRing)(n::n_Z) = n

(::SingularIntegerRing)(n::libSingular.number) = n_Z(n) 

function (R::SingularIntegerRing)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Void,
         (Ptr{BigInt}, Ptr{fmpz}), &a, &x)
   return R(libSingular.n_InitMPZ(a, R.ptr))   
end


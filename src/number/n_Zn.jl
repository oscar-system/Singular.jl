export ResidueRing, characteristic

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::SingularN_ZnRing) = n_Zn

parent(a::n_Zn) = a.parent

parent_type(::Type{n_Zn}) = SingularN_ZnRing

base_ring(a::n_Zn) = SingularZZ

base_ring(a::SingularN_ZnRing) = SingularZZ

function characteristic(R::SingularN_ZnRing)
   return SingularZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy(a::n_Zn)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::SingularN_ZnRing) = R(1)

zero(R::SingularN_ZnRing) = R(0)

function isone(n::n_Zn)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Zn)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

isunit(n::n_Zn) = gcd(n, parent(n)(characteristic(parent(n)))) == 1

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_Zn) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::SingularN_ZnRing)
   print(io, "Residue Ring of Integer Ring modulo ", characteristic(c))
end

function show(io::IO, n::n_Zn)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
   n.ptr = nn[]

   m = libSingular.StringEndS()
   s = unsafe_string(m) 
   libSingular.omFree(Ptr{Void}(m))

   print(io, s)
end

needs_parentheses(x::n_Zn) = false

isnegative(x::n_Zn) = false

show_minus_one(::Type{n_Zn}) = true

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Zn) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_Zn, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_Zn) = parent(y)(x) + y

-(x::n_Zn, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_Zn) = parent(y)(x) - y

*(x::n_Zn, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_Zn) = parent(y)(x)*y

+(x::n_Zn, y::n_Z) = x + parent(x)(y)

+(x::n_Z, y::n_Zn) = parent(y)(x) + y

-(x::n_Zn, y::n_Z) = x - parent(x)(y)

-(x::n_Z, y::n_Zn) = parent(y)(x) - y

*(x::n_Zn, y::n_Z) = x*parent(x)(y)

*(x::n_Z, y::n_Zn) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::n_Zn, y::n_Zn)
    libSingular.n_Greater(y.ptr, x.ptr, parent(x).ptr)
end

function ==(x::n_Zn, y::n_Zn)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end

isequal(x::n_Zn, y::n_Zn) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_Zn, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_Zn) = (parent(y)(x) == y)

==(x::n_Zn, y::n_Z) = (x ==  parent(x)(y))

==(x::n_Z, y::n_Zn) = (parent(y)(x) == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_Zn, y::Int)
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

function div(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_ExactDiv(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::n_Zn, y::n_Zn)
   par = parent(x)
   r = [libSingular.n_Init(0, par.ptr)]
   q = libSingular.n_QuotRem(x.ptr, y.ptr, pointer(r), par.ptr)
   return par(q), par(r)
end

function rem(x::n_Zn, y::n_Zn)
   par = parent(x)
   r = libSingular.n_IntMod(x.ptr, y.ptr, par.ptr)
   return par(r)
end

mod(x::n_Zn, y::n_Zn) = rem(x, y)

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_Zn, y::n_Zn)
   par = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

function lcm(x::n_Zn, y::n_Zn)
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

function gcdx(x::n_Zn, y::n_Zn)
   par = parent(x)
   s = [libSingular.n_Init(0, par.ptr)]
   t = [libSingular.n_Init(0, par.ptr)]
   g = libSingular.n_ExtGcd(x.ptr, y.ptr, pointer(s), pointer(t), par.ptr)
   return par(g), par(s[]), par(t[])
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Zn, y::n_Zn)
    xx = libSingular.number_ref(x.ptr)
    libSingular.n_InpAdd(xx, y.ptr, parent(x).ptr)
    x.ptr = xx[]
    nothing
end

function mul!(x::n_Zn, y::n_Zn, z::n_Zn)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function add!(x::n_Zn, y::n_Zn, z::n_Zn)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function zero!(x::n_Zn)
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

promote_rule{T <: Integer}(C::Type{n_Zn}, ::Type{T}) = n_Zn

promote_rule(C::Type{n_Zn}, ::Type{n_Z}) = n_Zn

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::SingularN_ZnRing)()
   z = n_Zn(R)
   z.parent = R
   return z
end

function (R::SingularN_ZnRing)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 
   z.parent = R
   return z
end

function (R::SingularN_ZnRing)(n::Int)
   z = n_Zn(R, n)
   z.parent = R
   return z
end

function (R::SingularN_ZnRing)(n::n_Z)
   m = libSingular.nApplyMapFunc(R.from_n_Z, n.ptr, parent(n).ptr, R.ptr)
   z = n_Zn(R, m)
   z.parent = R
   return z
end

(R::SingularN_ZnRing)(n::n_Zn) = n

function (R::SingularN_ZnRing)(n::libSingular.number)
   z = n_Zn(R, n) 
   z.parent = R
   return z
end

function (R::SingularN_ZnRing)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Void,
         (Ptr{BigInt}, Ptr{fmpz}), &a, &x)
   z = R(libSingular.n_InitMPZ(a, R.ptr))
      z.parent = R
   return z
end

###############################################################################
#
#   ResidueRing constructor
#
###############################################################################

function ResidueRing(R::SingularIntegerRing, a::Int; cached=true)
   a == 0 && throw(DivideError())
   a < 0 && throw(DomainError())

   return SingularN_ZnRing(a)
end

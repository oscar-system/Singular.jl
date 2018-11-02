export n_Z, Integers, crt

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{Integers}) = n_Z

parent(a::n_Z) = ZZ

parent_type(::Type{n_Z}) = Integers

base_ring(a::n_Z) = Union{}

base_ring(a::Integers) = Union{}

function deepcopy_internal(a::n_Z, dict::ObjectIdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(::Integers) = ZZ(1)

zero(::Integers) = ZZ(0)

function isone(n::n_Z)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Z)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
    isunit(n::n_Z)
> Return `true` if $n$ is $\pm 1$.
"""
isunit(n::n_Z) = n == 1 || n == -1

@doc Markdown.doc"""
    numerator(n::n_Z)
> Return the numerator of $n$ (which is $n$ itself).
"""
function numerator(n::n_Z)
   par = parent(n)
   r = libSingular.n_Copy(n.ptr, par.ptr)
   return par(r)
end

@doc Markdown.doc"""
    denominator(n::n_Z)
> Return the denominator of $n$ (which will always be $1$).
"""
function denominator(n::n_Z)
   return one(parent(n))
end

@doc Markdown.doc"""
    abs(n::n_Z)
> Return the absolute value of $n$.
"""
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

function show(io::IO, c::Integers)
   print(io, "Integer Ring")
end

function show(io::IO, n::n_Z)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
   n.ptr = nn

   m = libSingular.StringEndS()

   print(io,m)

end

needs_parentheses(x::n_Z) = false

isnegative(x::n_Z) = !libSingular.n_GreaterZero(x.ptr, parent(x).ptr) &&
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

function mod(x::n_Z, y::n_Z)
   d = rem(x, y)
   return d < 0 ? d + abs(y) : d
end

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
   println("type par",typeof(par.ptr))
   println("type r",typeof(r))
   println("type m",typeof(m))
   a = libSingular.n_ChineseRemainderSym(reinterpret(Ptr{Nothing},pointer(r)), reinterpret(Ptr{Nothing},pointer(m)), Cint(2), Cint(signed), par.ptr)
   println("Got through")
   return par(a)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Z, y::n_Z)
    x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
    return x
end

function mul!(x::n_Z, y::n_Z, z::n_Z)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_Z, y::n_Z, z::n_Z)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_Z)
   ptr = libSingular.n_Init(0, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule(C::Type{n_Z}, ::Type{T}) where {T <: Integer} = n_Z

###############################################################################
#
#   Parent call functions
#
###############################################################################

(::Integers)() = n_Z()

(R::Integers)(x::Integer) = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 

(::Integers)(n::Int) = n_Z(n)

(::Integers)(n::n_Z) = n

(::Integers)(n::libSingular.number) = n_Z(n) 

function (R::Integers)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Nothing,
         (Ptr{BigInt}, Ptr{fmpz}), &a, &x)
   return R(libSingular.n_InitMPZ(a, R.ptr))   
end


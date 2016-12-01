export reconstruct

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::SingularRationalField) = n_Q

parent(a::n_Q) = SingularQQ

parent_type(::Type{n_Q}) = SingularRationalField

base_ring(a::n_Q) = SingularZZ

base_ring(a::SingularRationalField) = SingularZZ

function deepcopy(a::n_Q)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(::SingularRationalField) = SingularQQ(1)

zero(::SingularRationalField) = SingularQQ(0)

function isone(n::n_Q)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Q)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

isunit(n::n_Q) = !iszero(n)

function num(n::n_Q)
   nn = libSingular.number_ref(n.ptr);
   p = libSingular.n_GetNumerator(nn, parent(n).ptr)
   pp = libSingular.nApplyMapFunc(n_Q_2_n_Z, p, SingularQQ.ptr, SingularZZ.ptr)
   libSingular.n_Delete(p, SingularQQ.ptr)
   return SingularZZ(pp)
end

function den(n::n_Q)
   nn = libSingular.number_ref(n.ptr);
   p = libSingular.n_GetDenom(nn, parent(n).ptr)
   pp = libSingular.nApplyMapFunc(n_Q_2_n_Z, p, SingularQQ.ptr, SingularZZ.ptr)
   libSingular.n_Delete(p, SingularQQ.ptr)
   return SingularZZ(pp)
end

function abs(n::n_Q)
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

canonical_unit(x::n_Q) = isnegative(x) ? -one(parent(x)) : one(parent(x))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::SingularRationalField)
   print(io, "Rational Field")
end

function show(io::IO, n::n_Q)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
   n.ptr = nn[]

   m = libSingular.StringEndS()
   s = bytestring(m) 
   libSingular.omFree(Ptr{Void}(m))

   print(io, s)
end

needs_parentheses(x::n_Q) = false

is_negative(x::n_Q) = !libSingular.n_GreaterZero(x.ptr, parent(x).ptr) &&
                      !iszero(x)

show_minus_one(::Type{n_Q}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Q) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_Q, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_Q) = parent(y)(x) + y

-(x::n_Q, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_Q) = parent(y)(x) - y

*(x::n_Q, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_Q) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::n_Q, y::n_Q)
    libSingular.n_Greater(y.ptr, x.ptr, parent(x).ptr)
end

function ==(x::n_Q, y::n_Q)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end

isequal(x::n_Q, y::n_Q) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_Q, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_Q) = (parent(y)(x) == y)

isequal(x::n_Q, y::Integer) = (x == y)

isequal(x::Integer, y::n_Q) = (x == y)

isless(x::n_Q, y::Integer) = isless(x, parent(x)(y))

isless(x::Integer, y::n_Q) = isless(parent(y)(x), y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_Q, y::Int)
    if y < 0
       return div(parent(x)(1), x^(-y))
    elseif isone(x)
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

function inv(x::n_Q)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function div(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::n_Q, y::n_Q)
   return div(x, y), zero(parent(x)) 
end

function rem(x::n_Q, y::n_Q)
   return zero(parent(x))
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_Q, y::n_Q)
   par = parent(x)
   p = libSingular.n_SubringGcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

function lcm(x::n_Q, y::n_Q)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   return div(x*y, gcd(x, y))
end

###############################################################################
#
#   Rational reconstruction
#
###############################################################################

function reconstruct(x::n_Z, y::n_Z)
   p = libSingular.n_Farey(x.ptr, y.ptr, SingularQQ.ptr)
   return SingularQQ(p)
end

reconstruct(x::n_Z, y::Integer) = reconstruct(x, SingularZZ(y))

reconstruct(x::Integer, y::n_Z) = reconstruct(SingularZZ(x), y)

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Q, y::n_Q)
    xx = libSingular.number_ref(x.ptr)
    libSingular.n_InpAdd(xx, y.ptr, parent(x).ptr)
    x.ptr = xx[]
    nothing
end

function mul!(x::n_Q, y::n_Q, z::n_Q)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function add!(x::n_Q, y::n_Q, z::n_Q)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function zero!(x::n_Q)
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

Base.promote_rule{T <: Integer}(C::Type{n_Q}, ::Type{T}) = n_Q

Base.promote_rule(C::Type{n_Q}, ::Type{Nemo.fmpz}) = n_Q

Base.promote_rule(C::Type{n_Q}, ::Type{n_Q}) = n_Z

###############################################################################
#
#   Parent call functions
#
###############################################################################

(::SingularRationalField)() = n_Q()

(::SingularRationalField)(n::Int) = n_Q(n)

(::SingularRationalField)(n::n_Z) = n_Q(n)

(::SingularRationalField)(n::libSingular.number) = n_Q(n) 

function (R::SingularRationalField)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Void,
         (Ptr{BigInt}, Ptr{fmpz}), &a, &x)
   return R(libSingular.n_InitMPZ(a, R.ptr))   
end


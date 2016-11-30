export gcd, lcm

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::SingularRationalField) = SingularQQElem

parent(a::SingularQQElem) = SingularQQ

parent_type(::Type{SingularQQElem}) = SingularRationalField

function check_parent(a::SingularQQElem, b::SingularQQElem) 
   parent(a) != parent(b) && error("Elements have different parents")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(::SingularRationalField) = SingularQQ(1)

zero(::SingularRationalField) = SingularQQ(0)

function isone(n::SingularQQElem)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::SingularQQElem)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::SingularRationalField)
   print(io, "Rational Field")
end

function show(io::IO, n::SingularQQElem)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
   n.ptr = nn[] ## TODO: FIXME: unsafe...

   m = libSingular.StringEndS()
   s = bytestring(m) 
   libSingular.omFree(Ptr{Void}(m))

   print(io, s)
end

needs_parentheses(x::SingularQQElem) = false

is_negative(x::SingularQQElem) = x < 0

show_minus_one(::Type{SingularQQElem}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::SingularQQElem) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::SingularQQElem, y::SingularQQElem)
   check_parent(x, y)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::SingularQQElem, y::SingularQQElem)
   check_parent(x, y)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::SingularQQElem, y::SingularQQElem)
   check_parent(x, y)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::SingularQQElem, y::Integer) = x + parent(x)(y)

+(x::Integer, y::SingularQQElem) = parent(y)(x) + y

-(x::SingularQQElem, y::Integer) = x - parent(x)(y)

-(x::Integer, y::SingularQQElem) = parent(y)(x) - y

*(x::SingularQQElem, y::Integer) = x*parent(x)(y)

*(x::Integer, y::SingularQQElem) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::SingularQQElem, y::SingularQQElem)
    check_parent(x, y)
    libSingular.n_Greater(y.ptr, x.ptr, parent(x).ptr)
end

function ==(x::SingularQQElem, y::SingularQQElem)
    check_parent(x, y)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end

isequal(x::SingularQQElem, y::SingularQQElem) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::SingularQQElem, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::SingularQQElem) = (parent(y)(x) == y)

isequal(x::SingularQQElem, y::Integer) = (x == y)

isequal(x::Integer, y::SingularQQElem) = (x == y)

isless(x::SingularQQElem, y::Integer) = isless(x, parent(x)(y))

isless(x::Integer, y::SingularQQElem) = isless(parent(y)(x), y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::SingularQQElem, y::Int)
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

function inv(x::SingularQQElem)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function div(x::SingularQQElem, y::SingularQQElem)
   check_parent(x, y)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function divexact(x::SingularQQElem, y::SingularQQElem)
   check_parent(x, y)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::SingularQQElem, y::SingularQQElem)
   check_parent(x, y)
   c = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function lcm(x::SingularQQElem, y::SingularQQElem)
   check_parent(x, y)
   c = parent(x)
   p = libSingular.n_Lcm(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::SingularQQElem, y::SingularQQElem)
    xx = libSingular.number_ref(x.ptr)
    libSingular.n_InpAdd(xx, y.ptr, parent(x).ptr)
    x.ptr = xx[]
    nothing
end

function mul!(x::SingularQQElem, y::SingularQQElem, z::SingularQQElem)
    ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
    libSingular.n_Delete(x.ptr, parent(x).ptr)
    x.ptr = ptr
    nothing
end

###############################################################################
#
#   Parent call functions
#
###############################################################################

(::SingularRationalField)() = SingularQQElem()

(::SingularRationalField)(n::Int) = SingularQQElem(n)

(::SingularRationalField)(n::libSingular.number) = SingularQQElem(n) 




using Nemo

###############################################################################
#
#   Memory management
#
###############################################################################

function fmpzInit(i::Clong, cf::Ptr{Cvoid})
   return number(Nemo.fmpz(i))
end
   
function fmpzDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
    ptr_new = unsafe_load(ptr)
    n = julia(ptr_new)::Nemo.fmpz
    if n != C_NULL
        number_pop!(nemoNumberID, ptr_new)
    end
    nothing
end

function fmpzCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fmpzGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    return Cint(n > 0)
end

function fmpzCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
    data_ptr = get_coeff_data(cf)
    r = unsafe_pointer_to_objref(data_ptr)
    str = string(r)
    libSingular.PrintS(str);
    nothing
end

function fmpzWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    if needs_parentheses(n)
        str = "("*string(n)*")"
    else
        str = string(n)
    end
    libSingular.StringAppendS(str);
    nothing
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function fmpzNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    return number(-n)
end

function fmpzInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    nptr = pointer_from_objref(n)
    ccall((:fmpz_neg, :libflint), Cvoid, (Ptr{Nemo.fmpz}, Ptr{Nemo.fmpz}), nptr, nptr)
    return number(n, false)
end

function fmpzInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    return number(divexact(1, n))
end

function fmpzMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n1 = julia(a)::Nemo.fmpz
    n2 = julia(b)::Nemo.fmpz
    return number(n1*n2)
end

function fmpzInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    r = unsafe_load(a)
    aa = julia(r)::Nemo.fmpz
    bb = julia(b)::Nemo.fmpz
    cc = mul!(aa, aa, bb)
    n = number(cc, aa)
    unsafe_store!(a, n, 1)
    nothing
end

function fmpzAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n1 = julia(a)::Nemo.fmpz
    n2 = julia(b)::Nemo.fmpz
    return number(n1 + n2)
end

function fmpzInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    r = unsafe_load(a)
    aa = julia(r)::Nemo.fmpz
    bb = julia(b)::Nemo.fmpz
    cc = addeq!(aa, bb)
    n = number(cc, aa)
    unsafe_store!(a, n, 1)
    nothing
end

function fmpzSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n1 = julia(a)::Nemo.fmpz
    n2 = julia(b)::Nemo.fmpz
    return number(n1 - n2)
end

function fmpzDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n1 = julia(a)::Nemo.fmpz
    n2 = julia(b)::Nemo.fmpz
    return number(divexact(n1, n2))
end

function fmpzDivBy(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n1 = julia(a)::Nemo.fmpz
    n2 = julia(b)::Nemo.fmpz
    return Cint(divides(n1, n2)[1])
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fmpzGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n1 = julia(a)::Nemo.fmpz
    n2 = julia(b)::Nemo.fmpz
    return Cint(n1 > n2)
end

function fmpzEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fmpz
   n2 = julia(b)::Nemo.fmpz
   return Cint(n1 == n2)
end

function fmpzIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fmpz
   return Cint(iszero(n))
end

function fmpzIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    return Cint(isone(n))
end

function fmpzIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n = julia(a)::Nemo.fmpz
    return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function fmpzGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
    n1 = julia(a)::Nemo.fmpz
    n2 = julia(b)::Nemo.fmpz
    return number(gcd(n1, n2))
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function fmpzExtGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, s::Ptr{Ptr{Cvoid}}, t::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fmpz
   n2 = julia(b)::Nemo.fmpz
   s1 = unsafe_load(s)
   if s1 != C_NULL
      number_pop!(nemoNumberID, s1)
   end
   t1 = unsafe_load(t)
   if t1 != C_NULL
      number_pop!(nemoNumberID, t1)
   end
   g1, s1, t1 = gcdx(n1, n2)
   setindex_internal_void(reinterpret(Ptr{Cvoid},s), number(s1))
   setindex_internal_void(reinterpret(Ptr{Cvoid},t), number(t1))
   return number(g1)
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fmpzInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
    n = julia(unsafe_load(ptr))::fmpz
    return Clong(n)
end

function fmpzMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
    n = julia(unsafe_load(ptr))::fmpz
    z = convert(BigInt, n)
    bptr = reinterpret(Ptr{Cvoid},pointer_from_objref(b))
    zptr = reinterpret(Ptr{Cvoid},pointer_from_objref(z))
    mpz_init_set_internal(bptr,zptr)
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function fmpzInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})

    ring_struct = singular_coeff_ring_struct()
    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 0
    ring_struct.is_domain = 1
    ring_struct.ch = 0
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(fmpzInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(fmpzInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(fmpzMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(fmpzInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(fmpzCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(fmpzDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(fmpzAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(fmpzInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(fmpzSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(fmpzMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(fmpzInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(fmpzDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDivBy = @cfunction(fmpzDivBy, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(fmpzInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(fmpzGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfExtGcd = @cfunction(fmpzExtGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(fmpzGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(fmpzEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(fmpzIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(fmpzIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(fmpzIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(fmpzGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(fmpzWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(fmpzCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct,cf)

    return Cint(0)
end

function register(R::FlintIntegerRing)
   c = @cfunction(fmpzInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   return nRegister(n_unknown, c)
end

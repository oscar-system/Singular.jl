###############################################################################
#
#   Memory management
#
###############################################################################

function fmpqInit(i::Clong, cf::Ptr{Cvoid})
   return number(Nemo.QQFieldElem(i))
end

function fmpqDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   ptr_new = unsafe_load(ptr)
   number_pop!(nemoNumberID, ptr_new)
   nothing
end

function fmpqCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fmpqGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   return Cint(n > 0)
end

function fmpqCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
    data_ptr = get_coeff_data(cf)
    r = unsafe_pointer_to_objref(data_ptr)
    str = string(r)
    libSingular.PrintS(str);
    nothing
end

function fmpqWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   libSingular.StringAppendS(AbstractAlgebra.obj_to_string_wrt_times(n))
   nothing
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function fmpqNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   return number(-n)
end

function fmpqInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
    ccall((:fmpq_neg, libflint), Cvoid, (Ptr{Nemo.QQFieldElem}, Ptr{Nemo.QQFieldElem}), a, a)
    return a
end

function fmpqInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   return number(Nemo.inv(n))
end

function fmpqMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.QQFieldElem
   n2 = julia(b)::Nemo.QQFieldElem
   return number(n1*n2)
end

function fmpqInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.QQFieldElem
   bb = julia(b)::Nemo.QQFieldElem
   Nemo.mul!(aa, aa, bb)
   nothing
end

function fmpqAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.QQFieldElem
   n2 = julia(b)::Nemo.QQFieldElem
   return number(n1 + n2)
end

function fmpqInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.QQFieldElem
   bb = julia(b)::Nemo.QQFieldElem
   Nemo.addeq!(aa, bb)
   nothing
end

function fmpqSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.QQFieldElem
   n2 = julia(b)::Nemo.QQFieldElem
   return number(n1 - n2)
end

function fmpqDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.QQFieldElem
   n2 = julia(b)::Nemo.QQFieldElem
   return number(Nemo.divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fmpqGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.QQFieldElem
   n2 = julia(b)::Nemo.QQFieldElem
   return Cint(n1 > n2)
end

function fmpqEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.QQFieldElem
   n2 = julia(b)::Nemo.QQFieldElem
   return Cint(n1 == n2)
end

function fmpqIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   return Cint(Nemo.iszero(n))
end

function fmpqIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   return Cint(Nemo.isone(n))
end

function fmpqIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.QQFieldElem
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function fmpqGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.QQFieldElem
   n2 = julia(b)::Nemo.QQFieldElem
   return number(Nemo.gcd(n1, n2))
end

###############################################################################
#
#   Subring GCD
#
###############################################################################

function fmpqSubringGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = numerator(julia(a)::Nemo.QQFieldElem)
   n2 = numerator(julia(b)::Nemo.QQFieldElem)
   return number(Nemo.QQFieldElem(Nemo.gcd(n1, n2)))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fmpqInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
    ptr_load = unsafe_load(ptr)
    n = julia(ptr_load)::Nemo.QQFieldElem
    ret_val = Clong(n)
    number_pop!(nemoNumberID, ptr_load)
    return ret_val
end

function fmpqMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
    ptr_load = unsafe_load(ptr)
    n = julia(unsafe_load(ptr))::Nemo.ZZRingElem
    z = convert(BigInt, Nemo.numerator(n))
    GC.@preserve b z begin
        bptr = reinterpret(Ptr{Cvoid}, pointer_from_objref(b))
        zptr = reinterpret(Ptr{Cvoid}, pointer_from_objref(z))
        number_pop!(nemoNumberID, ptr_load)
        mpz_init_set_internal(bptr, zptr)
    end
    nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function fmpqInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})
    ring_struct = singular_coeff_ring_struct()
    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 1
    ring_struct.is_domain = 1
    ring_struct.ch = 0
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(fmpqInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(fmpqInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(fmpqMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(fmpqInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(fmpqCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(fmpqDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(fmpqAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(fmpqInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(fmpqSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(fmpqMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(fmpqInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(fmpqDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(fmpqInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(fmpqGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSubringGcd = @cfunction(fmpqSubringGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(fmpqGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(fmpqEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(fmpqIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(fmpqIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(fmpqIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(fmpqGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(fmpqWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(fmpqCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct, cf)

    return Cint(0)
end

function register(R::Nemo.QQField)
   c = @cfunction(fmpqInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   return nRegister(n_unknown, c)
end

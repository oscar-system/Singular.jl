###############################################################################
#
#   Memory management
#
###############################################################################

function fqInit(i::Clong, cf::Ptr{Cvoid})
   cf_ptr = get_coeff_data_void(cf)
   R = julia(cf_ptr)::Nemo.FqFiniteField
   return number(R(Int(i)))
end

function fqDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   ptr_new = unsafe_load(ptr)
   number_pop!(nemoNumberID, Ptr{Cvoid}(ptr_new))
   nothing
end

function fqCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fqGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   return Cint(1)
end

function fqCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
   data_ptr = get_coeff_data(cf)
   r = unsafe_pointer_to_objref(data_ptr)
   str = string(r)
   libSingular.PrintS(str);
   nothing
end

function fqWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   libSingular.StringAppendS(libSingular.stringify_wrt_times(n))
   nothing
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function fqNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return number(-n)
end

function fqInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   cf_ptr = get_coeff_data_void(cf)
   ccall((:fq_neg, libflint), Cvoid,
	 (Ptr{Nemo.fq}, Ptr{Nemo.fq}, Ptr{Nemo.FqFiniteField}),
	   a, a, cf_ptr)
   return a
end

function fqInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return number(Nemo.inv(n))
end

function fqMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1*n2)
end

function fqInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq
   bb = julia(b)::Nemo.fq
   Nemo.mul!(aa, aa, bb)
   nothing
end

function fqAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1 + n2)
end

function fqInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq
   bb = julia(b)::Nemo.fq
   Nemo.addeq!(aa, bb)
   nothing
end

function fqSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1 - n2)
end

function fqDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(Nemo.divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fqGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return Cint(n1 != n2)
end

function fqEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return Cint(n1 == n2)
end

function fqIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return Cint(Nemo.iszero(n))
end

function fqIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return Cint(Nemo.isone(n))
end

function fqIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function fqGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(Nemo.gcd(n1, n2))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fqInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   return Clong(0)
end

function fqMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   ptr_load = unsafe_load(ptr)
   z = BigInt(0)
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

function fqInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})

    R = julia(p)

    ring_struct = singular_coeff_ring_struct()
    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 1
    ring_struct.is_domain = 1
    ring_struct.ch = Cint(BigInt(Nemo.characteristic(R)))
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(fqInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(fqInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(fqMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(fqInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(fqCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(fqDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(fqAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(fqInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(fqSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(fqMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(fqInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(fqDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(fqInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(fqGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(fqGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(fqEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(fqIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(fqIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(fqIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(fqGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(fqWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(fqCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct, cf)

    return Cint(0)
end

function register(R::Nemo.FqFiniteField)
   c = @cfunction(fqInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   return nRegister(n_unknown, c)
end

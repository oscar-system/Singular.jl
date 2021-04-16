###############################################################################
#
#   Memory management
#
###############################################################################

function fq_nmodInit(i::Clong, cf::Ptr{Cvoid})
   cf_ptr = get_coeff_data_void(cf)
   R = julia(cf_ptr)::Nemo.FqNmodFiniteField
   return number(R(Int(i)))
end

function fq_nmodDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   ptr_new = unsafe_load(ptr)
   number_pop!(nemoNumberID, Ptr{Cvoid}(ptr_new))
   nothing
end

function fq_nmodCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fq_nmodGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   return Cint(1)
end

function fq_nmodCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
   data_ptr = get_coeff_data(cf)
   r = unsafe_pointer_to_objref(data_ptr)
   str = string(r)
   libSingular.PrintS(str);
   nothing
end

function fq_nmodWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   libSingular.StringAppendS(libSingular.stringify_wrt_times(n))
   nothing
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function fq_nmodNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return number(-n)
end

function fq_nmodInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   cf_ptr = get_coeff_data_void(cf)
   ccall((:fq_nmod_neg, libflint), Cvoid,
	 (Ptr{Nemo.fq_nmod}, Ptr{Nemo.fq_nmod}, Ptr{Nemo.FqNmodFiniteField}),
	   a, a, cf_ptr)
   return a
end

function fq_nmodInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return number(Nemo.inv(n))
end

function fq_nmodMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(n1*n2)
end

function fq_nmodInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq_nmod
   bb = julia(b)::Nemo.fq_nmod
   Nemo.mul!(aa, aa, bb)
   nothing
end

function fq_nmodAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(n1 + n2)
end

function fq_nmodInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq_nmod
   bb = julia(b)::Nemo.fq_nmod
   Nemo.addeq!(aa, bb)
   nothing
end

function fq_nmodSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(n1 - n2)
end

function fq_nmodDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(Nemo.divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fq_nmodGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return Cint(n1 != n2)
end

function fq_nmodEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return Cint(n1 == n2)
end

function fq_nmodIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return Cint(Nemo.iszero(n))
end

function fq_nmodIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return Cint(Nemo.isone(n))
end

function fq_nmodIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function fq_nmodGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(Nemo.gcd(n1, n2))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fq_nmodInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   return Clong(0)
end

function fq_nmodMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
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

function fq_nmodInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})

    R = julia(p)

    ring_struct = singular_coeff_ring_struct()
    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 1
    ring_struct.is_domain = 1
    ring_struct.ch = Cint(BigInt(Nemo.characteristic(R)))
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(fq_nmodInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(fq_nmodInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(fq_nmodMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(fq_nmodInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(fq_nmodCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(fq_nmodDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(fq_nmodAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(fq_nmodInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(fq_nmodSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(fq_nmodMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(fq_nmodInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(fq_nmodDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(fq_nmodInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(fq_nmodGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(fq_nmodGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(fq_nmodEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(fq_nmodIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(fq_nmodIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(fq_nmodIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(fq_nmodGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(fq_nmodWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(fq_nmodCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct, cf)

    return Cint(0)
end

function register(R::Nemo.FqNmodFiniteField)
   c = @cfunction(fq_nmodInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   return nRegister(n_unknown, c)
end

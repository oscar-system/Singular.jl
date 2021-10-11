###############################################################################
#
#   Memory management
#
###############################################################################

function nf_elemInit(i::Clong, cf::Ptr{Cvoid})
   cf_ptr = get_coeff_data_void(cf)
   R = julia(cf_ptr)::Nemo.AnticNumberField
   return number(R(Int(i)))
end

function nf_elemDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   ptr_new = unsafe_load(ptr)
   number_pop!(nemoNumberID, Ptr{Cvoid}(ptr_new))
   nothing
end

function nf_elemCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.nf_elem
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function nf_elemGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   return Cint(1)
end

function nf_elemCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
   data_ptr = get_coeff_data(cf)
   r = unsafe_pointer_to_objref(data_ptr)
   str = string(r)
   libSingular.PrintS(str);
   nothing
end

function nf_elemWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.nf_elem
   libSingular.StringAppendS(AbstractAlgebra.obj_to_string_wrt_times(n))
   nothing
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function nf_elemNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.nf_elem
   return number(-n)
end

function nf_elemInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   cf_ptr = get_coeff_data_void(cf)
   ccall((:nf_elem_neg, libantic), Cvoid,
	 (Ptr{Nemo.nf_elem}, Ptr{Nemo.nf_elem}, Ptr{Nemo.AnticNumberField}),
	   a, a, cf_ptr)
   return a
end

function nf_elemInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.nf_elem
   return number(Nemo.inv(n))
end

function nf_elemMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(n1*n2)
end

function nf_elemInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.nf_elem
   bb = julia(b)::Nemo.nf_elem
   Nemo.mul!(aa, aa, bb)
   nothing
end

function nf_elemAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(n1 + n2)
end

function nf_elemInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.nf_elem
   bb = julia(b)::Nemo.nf_elem
   Nemo.addeq!(aa, bb)
   nothing
end

function nf_elemSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(n1 - n2)
end

function nf_elemDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(Nemo.divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function nf_elemGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return Cint(n1 != n2)
end

function nf_elemEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return Cint(n1 == n2)
end

function nf_elemIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.nf_elem
   return Cint(Nemo.iszero(n))
end

function nf_elemIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.nf_elem
   return Cint(Nemo.isone(n))
end

function nf_elemIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.nf_elem
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function nf_elemGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(Nemo.gcd(n1, n2))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function nf_elemInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   return Clong(0)
end

function nf_elemMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
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

function nf_elemInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})

    R = julia(p)

    ring_struct = singular_coeff_ring_struct()
    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 1
    ring_struct.is_domain = 1
    ring_struct.ch = 0
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(nf_elemInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(nf_elemInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(nf_elemMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(nf_elemInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(nf_elemCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(nf_elemDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(nf_elemAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(nf_elemInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(nf_elemSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(nf_elemMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(nf_elemInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(nf_elemDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(nf_elemInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(nf_elemGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(nf_elemGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(nf_elemEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(nf_elemIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(nf_elemIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(nf_elemIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(nf_elemGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(nf_elemWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(nf_elemCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct, cf)

    return Cint(0)
end

function register(R::Nemo.AnticNumberField)
   c = @cfunction(nf_elemInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   return nRegister(n_unknown, c)
end

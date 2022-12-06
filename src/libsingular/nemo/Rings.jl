###############################################################################
#
#   Memory management
#
###############################################################################

function nemoRingInit(i::Clong, cf::Ptr{Cvoid})
   data_ptr = get_coeff_data_void(cf)
   R = unsafe_pointer_to_objref(data_ptr)
   return number(R(i))
end

function nemoRingDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, n)
   end
   nothing
end

function nemoRingCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function nemoRingGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   ## Fixme: Is this in any sense correct?
   return Cint(1)
end

function nemoRingCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
  data_ptr = get_coeff_data(cf)
  r = unsafe_pointer_to_objref(data_ptr)
  str = string(r)
  libSingular.PrintS(str);
  nothing
end

function nemoRingWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   libSingular.StringAppendS(AbstractAlgebra.obj_to_string_wrt_times(n))
   nothing
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function nemoRingNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(-n)
end

function nemoRingInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   data_ptr = get_coeff_data_void(cf)
   R = unsafe_pointer_to_objref(data_ptr)
   n = julia(a)
   mone = R(-1)
   n_new = Nemo.mul!(n, n, mone)
   return number(n_new, n)
end

function nemoRingInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(Nemo.divexact(1, n))
end

function nemoRingMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1*n2)
end

function nemoRingInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   cc = Nemo.mul!(aa, aa, bb)
   n = number(cc, aa)
   setindex_internal_void(reinterpret(Ptr{Cvoid},a), n)
   nothing
end

function nemoRingAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 + n2)
end

function nemoRingInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   cc = Nemo.addeq!(aa, bb)
   n = number(cc, aa)
   setindex_internal_void(reinterpret(Ptr{Cvoid},a), n)
   nothing
end

function nemoRingSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 - n2)
end

function nemoRingDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(Nemo.divexact(n1, n2))
end

function nemoRingDivBy(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   if Nemo.iszero(n1)
     return Cint(Nemo.is_zero_divisor(n2))
   end
   return Cint(Nemo.divides(n1, n2)[1])
end

###############################################################################
#
#   Comparison
#
###############################################################################

function nemoRingGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 != n2)
end

function nemoRingEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 == n2)
end

function nemoRingIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(Nemo.iszero(n))
end

function nemoRingIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(Nemo.isone(n))
end

function nemoRingIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(n == -1)
end

###############################################################################
#
#   Annihilator
#
###############################################################################

function nemoRingAnn(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   f, b = AbstractAlgebra.is_zero_divisor_with_annihilator(n)
   return f ? number(b) : number(zero(parent(n)))
end

function nemoDomainAnn(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(iszero(n) ? one(parent(n)) : zero(parent(n)))
end

###############################################################################
#
#   GCD
#
###############################################################################

function nemoRingGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(Nemo.gcd(n1, n2))
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function nemoRingExtGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, s::Ptr{Ptr{Cvoid}}, t::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.RingElem
   n2 = julia(b)::Nemo.RingElem
   g1, s1, t1 = Nemo.gcdx(n1, n2)
   setindex_internal_void(reinterpret(Ptr{Cvoid}, s), number(s1))
   setindex_internal_void(reinterpret(Ptr{Cvoid}, t), number(t1))
   return number(g1)
end

###############################################################################
#
#   Conversion
#
###############################################################################

function nemoRingInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   return Clong(0)
end

function nemoRingMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   GC.@preserve b begin
      bptr = pointer_from_objref(b)
      mpz_init_set_si_internal(reinterpret(Ptr{Cvoid}, bptr), 0)
   end
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################
function nemoRingInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})

    ring_struct = singular_coeff_ring_struct()

    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 0
    ring_struct.is_domain = 0
    ring_struct.ch = 0
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(nemoRingInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(nemoRingInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(nemoRingMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(nemoRingInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(nemoRingCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(nemoRingDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(nemoRingAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(nemoRingInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(nemoRingSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(nemoRingMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(nemoRingInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(nemoRingDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDivBy = @cfunction(nemoRingDivBy, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(nemoRingInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfAnn = @cfunction(nemoRingAnn, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(nemoRingGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfExtGcd = @cfunction(nemoRingExtGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(nemoRingGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(nemoRingEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(nemoRingIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(nemoRingIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(nemoRingIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(nemoRingGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(nemoRingWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(nemoRingCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct,cf)

    return Cint(0)
end

function nemoDomainInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})

    ring_struct = singular_coeff_ring_struct()

    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 0
    ring_struct.is_domain = 1
    ring_struct.ch = 0
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(nemoRingInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(nemoRingInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(nemoRingMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(nemoRingInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(nemoRingCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(nemoRingDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(nemoRingAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(nemoRingInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(nemoRingSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(nemoRingMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(nemoRingInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(nemoRingDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDivBy = @cfunction(nemoRingDivBy, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(nemoRingInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfAnn = @cfunction(nemoDomainAnn, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(nemoRingGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfExtGcd = @cfunction(nemoRingExtGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(nemoRingGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(nemoRingEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(nemoRingIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(nemoRingIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(nemoRingIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(nemoRingGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(nemoRingWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(nemoRingCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct,cf)

    return Cint(0)
end

function register(R::Nemo.Ring)
   if AbstractAlgebra.is_domain_type(AbstractAlgebra.elem_type(R))
      c = @cfunction(nemoDomainInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   else
      c = @cfunction(nemoRingInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   end
   return nRegister(n_unknown, c)
end

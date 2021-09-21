###############################################################################
#
#   Memory management
#
###############################################################################

function nemoFieldInit(i::Clong, cf::Ptr{Cvoid})
   data_ptr = get_coeff_data_void(cf)
   R = unsafe_pointer_to_objref(data_ptr)
   return number(R(i))
end

function nemoFieldDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, n)
   end
   nothing
end

function nemoFieldCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function nemoFieldGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   ## Fixme: Is this in any sense correct?
   return Cint(1)
end

function nemoFieldCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
  data_ptr = get_coeff_data(cf)
  r = unsafe_pointer_to_objref(data_ptr)
  str = string(r)
  libSingular.PrintS(str);
  nothing
end

function nemoFieldWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   libSingular.StringAppendS(AbstractAlgebra.obj_to_string_wrt_times(n))
   nothing
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function nemoFieldNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(-n)
end

function nemoFieldInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   data_ptr = get_coeff_data_void(cf)
   R = unsafe_pointer_to_objref(data_ptr)
   n = julia(a)
   mone = R(-1)
   n_new = Nemo.mul!(n, n, mone)
   return number(n_new, n)
end

function nemoFieldInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(Nemo.inv(n))
end

function nemoFieldMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1*n2)
end

function nemoFieldInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   cc = Nemo.mul!(aa, aa, bb)
   n = number(cc, aa)
   setindex_internal_void(reinterpret(Ptr{Cvoid},a), n)
   nothing
end

function nemoFieldAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 + n2)
end

function nemoFieldInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   cc = Nemo.addeq!(aa, bb)
   n = number(cc, aa)
   setindex_internal_void(reinterpret(Ptr{Cvoid},a), n)
   nothing
end

function nemoFieldSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 - n2)
end

function nemoFieldDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(Nemo.divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function nemoFieldGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 != n2)
end

function nemoFieldEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 == n2)
end

function nemoFieldIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(Nemo.iszero(n))
end

function nemoFieldIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(Nemo.isone(n))
end

function nemoFieldIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function nemoFieldGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(Nemo.gcd(n1, n2))
end

###############################################################################
#
#   Subring GCD
#
###############################################################################

function nemoFieldSubringGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   data_ptr = get_coeff_data_void(cf)
   R =unsafe_load(data_ptr)
   if isa(R, Nemo.FracField)
      n1 = numerator(julia(a))
      n2 = numerator(julia(b))
      return number(R(Nemo.gcd(n1, n2)))
   else
      return number(deepcopy(julia(a)))
   end
end

###############################################################################
#
#   Conversion
#
###############################################################################

function nemoFieldInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   return Clong(0)
end

function nemoFieldMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   GC.@preserve b begin
      bptr = pointer_from_objref(b)
      mpz_init_set_si_internal(reinterpret(Ptr{Cvoid}, bptr),0)
   end
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function nemoFieldInitChar(cf::Ptr{Cvoid}, p::Ptr{Cvoid})

    ring_struct = singular_coeff_ring_struct()

    ring_struct.has_simple_alloc = 0
    ring_struct.has_simple_inverse = 0
    ring_struct.is_field = 1
    ring_struct.is_domain = 1
    ring_struct.ch = 0
    ring_struct.data = p
    ring_struct.cfInit = @cfunction(nemoFieldInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    ring_struct.cfInt = @cfunction(nemoFieldInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfMPZ = @cfunction(nemoFieldMPZ, Cvoid, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfInpNeg = @cfunction(nemoFieldInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCopy = @cfunction(nemoFieldCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDelete = @cfunction(nemoFieldDelete, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    ring_struct.cfAdd = @cfunction(nemoFieldAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpAdd = @cfunction(nemoFieldInpAdd, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSub = @cfunction(nemoFieldSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfMult = @cfunction(nemoFieldMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInpMult = @cfunction(nemoFieldInpMult, Cvoid, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfDiv = @cfunction(nemoFieldDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfInvers = @cfunction(nemoFieldInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGcd = @cfunction(nemoFieldGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfSubringGcd = @cfunction(nemoFieldSubringGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreater = @cfunction(nemoFieldGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfEqual = @cfunction(nemoFieldEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsZero = @cfunction(nemoFieldIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsOne = @cfunction(nemoFieldIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfIsMOne = @cfunction(nemoFieldIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfGreaterZero = @cfunction(nemoFieldGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfWriteLong = @cfunction(nemoFieldWrite, Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}))
    ring_struct.cfCoeffWrite = @cfunction(nemoFieldCoeffWrite, Cvoid, (Ptr{Cvoid}, Cint))

    fill_coeffs_with_function_data(ring_struct,cf)

    return Cint(0)
end

function register(R::Nemo.Field)
   c = @cfunction(nemoFieldInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   return nRegister(n_unknown, c)
end

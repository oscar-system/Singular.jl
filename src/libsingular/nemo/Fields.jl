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
   cc = Nemo.add!(aa, bb)
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
#   Factorization
#
###############################################################################

function _unwrap_nemo_field_elem(a)
   if isdefined(Singular, :FieldElemWrapper) && a isa Singular.FieldElemWrapper
      return a.data
   end
   return a
end

function _nemo_field_for_factorization(R)
   K = Singular.base_ring(R).base_ring
   if isdefined(Singular, :FieldWrapper) && K isa Singular.FieldWrapper
      return K.data
   end
   return K
end

function _singular_univariate_to_julia(R, f)
   K = _nemo_field_for_factorization(R)
   U, x = AbstractAlgebra.polynomial_ring(K, string(R.S[1]))
   g = zero(U)
   for (c, v) in zip(Singular.coefficients(f), Singular.exponent_vectors(f))
      cc = GC.@preserve c julia(cast_number_to_void(c.ptr))
      g += _unwrap_nemo_field_elem(cc)*x^v[1]
   end
   return g
end

function _julia_univariate_to_singular(R, f)
   S = Singular.base_ring(R)
   B = Singular.MPolyBuildCtx(R)
   for i in 0:Nemo.degree(f)
      c = Nemo.coeff(f, i)
      if !Nemo.iszero(c)
         Singular.push_term!(B, S(c), [i])
      end
   end
   return Singular.finish(B)
end

function _factorization_callback_ring(r_ptr::Ptr{Cvoid})
   r = factorization_callback_ring_copy(r_ptr)
   K = unsafe_pointer_to_objref(factorization_callback_coeff_data(r_ptr))
   S = Singular.CoefficientRing(K)
   T = Singular.elem_type(S)
   return Singular.PolyRing{T}(r, S, Singular.singular_symbols(r))
end

function nemoFieldFactorize(f_ptr::Ptr{Cvoid}, v_ptr::Ptr{Ptr{Cvoid}},
                            with_exps::Cint, r_ptr::Ptr{Cvoid})
   try
      R = _factorization_callback_ring(r_ptr)
      if Singular.nvars(R) != 1
         return C_NULL
      end

      f = R(factorization_callback_poly_copy_to_ring(f_ptr, r_ptr, R.ptr))
      F = AbstractAlgebra.factor(_singular_univariate_to_julia(R, f))

      factors = typeof(f)[]
      exponents = Int32[]
      if with_exps == 0 || with_exps == 3
         push!(factors, _julia_univariate_to_singular(R, AbstractAlgebra.unit(F)))
         push!(exponents, Int32(1))
      end
      for (g, e) in F
         push!(factors, _julia_univariate_to_singular(R, g))
         push!(exponents, Int32(e))
      end

      I = Singular.Ideal(R, factors)
      res = GC.@preserve R f I copy_factorization_result(I.ptr, exponents,
         reinterpret(Ptr{Cvoid}, v_ptr), Int(with_exps), R.ptr, r_ptr)
      return reinterpret(Ptr{Cvoid}, res.cpp_object)
   catch
      return C_NULL
   end
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
    ring_struct.cfAnn = @cfunction(nemoDomainAnn, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
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
    ring_struct.cfFactorize = @cfunction(nemoFieldFactorize, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Cint, Ptr{Cvoid}))

    fill_coeffs_with_function_data(ring_struct,cf)

    return Cint(0)
end

function register(R::Nemo.Field)
   c = @cfunction(nemoFieldInitChar, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
   return nRegister(n_Nemo_Field, c)
end

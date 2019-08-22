###############################################################################
#
#   Memory management
#
###############################################################################

function fq_nmodInit(i::Clong, cf::Ptr{Cvoid})
   R = julia(cf)::Nemo.FqNmodFiniteField
   return number(R(Int(i)))
end
   
function fq_nmodDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
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
   r = julia(cf)::Nemo.FqNmodFiniteField
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function fq_nmodWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   if needs_parentheses(n)
      str = "("*string(n)*")"
   else
      str = string(n)
   end
   icxx"""StringAppendS($str);"""
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
   R = julia(cf)::Nemo.FqNmodFiniteField
   n = julia(a)::Nemo.fq_nmod
   ccall((:fq_nmod_neg, :libflint), Void, (Ptr{Nemo.fq_nmod}, Ptr{Nemo.fq_nmod}, Ptr{FqNmodFiniteField}), &n, &n, &R)
   return number(n, false)
end

function fq_nmodInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return number(inv(n))
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
   cc = mul!(aa, aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
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
   cc = addeq!(aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
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
   return number(divexact(n1, n2))
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
   return Cint(iszero(n))
end

function fq_nmodIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq_nmod
   return Cint(isone(n))
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
   return number(gcd(n1, n2))
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
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function fq_nmodInitChar(cf::Ptr{Cvoid}, p::Ptr{Void})
        
    pInit = cfunction(fq_nmodInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    pInt = cfunction(fq_nmodInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pMPZ = cfunction(fq_nmodMPZ, Void, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pInpNeg = cfunction(fq_nmodInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pCopy = cfunction(fq_nmodCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pDelete = cfunction(fq_nmodDelete, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pAdd = cfunction(fq_nmodAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInpAdd = cfunction(fq_nmodInpAdd, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    pSub = cfunction(fq_nmodSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pMult = cfunction(fq_nmodMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInpMult = cfunction(fq_nmodInpMult, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    pDiv = cfunction(fq_nmodDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInvers = cfunction(fq_nmodInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pGcd = cfunction(fq_nmodGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pGreater = cfunction(fq_nmodGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pEqual = cfunction(fq_nmodEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pIsZero = cfunction(fq_nmodIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pIsOne = cfunction(fq_nmodIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pIsMOne = cfunction(fq_nmodIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pGreaterZero = cfunction(fq_nmodGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pWrite = cfunction(fq_nmodWrite, Void, (Ptr{Cvoid}, Ptr{Cvoid}))
    pCoeffWrite = cfunction(fq_nmodCoeffWrite, Void, (Ptr{Cvoid}, Cint))

    return Cint(0)
end

function register(R::FqNmodFiniteField)
   c = cfunction(fq_nmodInitChar, Cint, (Ptr{Cvoid}, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end

###############################################################################
#
#   Memory management
#
###############################################################################

function fq_nmodInit(i::Clong, cf::coeffs)
   R = julia(cf)::Nemo.FqNmodFiniteField
   return number(R(Int(i)))
end
   
function fq_nmodDelete(ptr::Ptr{number}, cf::coeffs)
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function fq_nmodCopy(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq_nmod
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fq_nmodGreaterZero(a::number, cf::coeffs)
   return Cint(1)
end

function fq_nmodCoeffWrite(cf::coeffs, d::Cint)
   r = julia(cf)::Nemo.FqNmodFiniteField
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function fq_nmodWrite(a::number, cf::coeffs)
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

function fq_nmodNeg(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq_nmod
   return number(-n)
end

function fq_nmodInpNeg(a::number, cf::coeffs)
   R = julia(cf)::Nemo.FqNmodFiniteField
   n = julia(a)::Nemo.fq_nmod
   ccall((:fq_nmod_neg, :libflint), Void, (Ptr{Nemo.fq_nmod}, Ptr{Nemo.fq_nmod}, Ptr{FqNmodFiniteField}), &n, &n, &R)
   return number(n, false)
end

function fq_nmodInvers(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq_nmod
   return number(inv(n))
end

function fq_nmodMult(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(n1*n2)
end

function fq_nmodInpMult(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq_nmod
   bb = julia(b)::Nemo.fq_nmod
   ptr1 = pointer_from_objref(aa)
   aa = mul!(aa, aa, bb)
   ptr2 = pointer_from_objref(aa)
   n = number(aa, ptr1 != ptr2)
   unsafe_store!(a, n, 1)
   nothing
end

function fq_nmodAdd(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(n1 + n2)
end

function fq_nmodInpAdd(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq_nmod
   bb = julia(b)::Nemo.fq_nmod
   ptr1 = pointer_from_objref(aa)
   aa = addeq!(aa, bb)
   ptr2 = pointer_from_objref(aa)
   n = number(aa, ptr1 != ptr2)
   unsafe_store!(a, n, 1)
   nothing
end

function fq_nmodSub(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(n1 - n2)
end

function fq_nmodDiv(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fq_nmodGreater(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return Cint(n1 != n2)
end

function fq_nmodEqual(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return Cint(n1 == n2)
end

function fq_nmodIsZero(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq_nmod
   return Cint(iszero(n))
end

function fq_nmodIsOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq_nmod
   return Cint(isone(n))
end

function fq_nmodIsMOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq_nmod
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function fq_nmodGcd(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq_nmod
   n2 = julia(b)::Nemo.fq_nmod
   return number(gcd(n1, n2))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fq_nmodInt(ptr::Ptr{number}, cf::coeffs)
   return Clong(0)
end

function fq_nmodMPZ(b::BigInt, ptr::Ptr{number}, cf::coeffs)
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function fq_nmodInitChar(cf::coeffs, p::Ptr{Void})
        
    pInit = cfunction(fq_nmodInit, number, (Clong, coeffs))
    pInt = cfunction(fq_nmodInt, Clong, (Ptr{number}, coeffs))
    pMPZ = cfunction(fq_nmodMPZ, Void, (BigInt, Ptr{number}, coeffs))
    pInpNeg = cfunction(fq_nmodInpNeg, number, (number, coeffs))
    pCopy = cfunction(fq_nmodCopy, number, (number, coeffs))
    pDelete = cfunction(fq_nmodDelete, Void, (Ptr{number}, coeffs))
    pAdd = cfunction(fq_nmodAdd, number, (number, number, coeffs))
    pInpAdd = cfunction(fq_nmodInpAdd, Void, (Ptr{number}, number, coeffs))
    pSub = cfunction(fq_nmodSub, number, (number, number, coeffs))
    pMult = cfunction(fq_nmodMult, number, (number, number, coeffs))
    pInpMult = cfunction(fq_nmodInpMult, Void, (Ptr{number}, number, coeffs))
    pDiv = cfunction(fq_nmodDiv, number, (number, number, coeffs))
    pInvers = cfunction(fq_nmodInvers, number, (number, coeffs))
    pGcd = cfunction(fq_nmodGcd, number, (number, number, coeffs))
    pGreater = cfunction(fq_nmodGreater, Cint, (number, number, coeffs))
    pEqual = cfunction(fq_nmodEqual, Cint, (number, number, coeffs))
    pIsZero = cfunction(fq_nmodIsZero, Cint, (number, coeffs))
    pIsOne = cfunction(fq_nmodIsOne, Cint, (number, coeffs))
    pIsMOne = cfunction(fq_nmodIsMOne, Cint, (number, coeffs))
    pGreaterZero = cfunction(fq_nmodGreaterZero, Cint, (number, coeffs))
    pWrite = cfunction(fq_nmodWrite, Void, (number, coeffs))
    pCoeffWrite = cfunction(fq_nmodCoeffWrite, Void, (coeffs, Cint))

    icxx""" 
      coeffs cf = (coeffs)($cf);
      cf->has_simple_Alloc = FALSE;  
      cf->has_simple_Inverse= FALSE;          
      cf->is_field  = TRUE;
      cf->is_domain = TRUE;
      cf->ch = 0;
      cf->data = $p;
      cf->cfInit = (number (*)(long, const coeffs)) $pInit;
      cf->cfInt = (long (*)(number &, const coeffs)) $pInt;
      cf->cfMPZ = (void (*)(__mpz_struct *, number &, const coeffs)) $pMPZ;
      cf->cfInpNeg = (number (*)(number, const coeffs)) $pInpNeg;
      cf->cfCopy = (number (*)(number, const coeffs)) $pCopy;
      cf->cfDelete = (void (*)(number *, const coeffs)) $pDelete;
      cf->cfAdd = (numberfunc) $pAdd;
      cf->cfInpAdd = (void (*)(number &, number, const coeffs)) $pInpAdd;
      cf->cfSub = (numberfunc) $pSub;
      cf->cfMult = (numberfunc) $pMult;
      cf->cfInpMult = (void (*)(number &, number, const coeffs)) $pInpMult;
      cf->cfDiv = (numberfunc) $pDiv;
      cf->cfInvers = (number (*)(number, const coeffs)) $pInvers;
      cf->cfGcd = (numberfunc) $pGcd;
      cf->cfGreater = (BOOLEAN (*)(number, number, const coeffs)) $pGreater;
      cf->cfEqual = (BOOLEAN (*)(number, number, const coeffs)) $pEqual;
      cf->cfIsZero = (BOOLEAN (*)(number, const coeffs)) $pIsZero;
      cf->cfIsOne = (BOOLEAN (*)(number, const coeffs)) $pIsOne;
      cf->cfIsMOne = (BOOLEAN (*)(number, const coeffs)) $pIsMOne;
      cf->cfGreaterZero = (BOOLEAN (*)(number, const coeffs)) $pGreaterZero;
      cf->cfWriteLong = (void (*)(number, const coeffs)) $pWrite;
      cf->cfCoeffWrite = (void (*)(const coeffs, BOOLEAN)) $pCoeffWrite;
    """

    return Cint(0)
end

function register(R::FqNmodFiniteField)
   c = cfunction(fq_nmodInitChar, Cint, (coeffs, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end

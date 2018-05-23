###############################################################################
#
#   Memory management
#
###############################################################################

function nf_elemInit(i::Clong, cf::coeffs)
   R = julia(cf)::Nemo.AnticNumberField
   return number(R(Int(i)))
end
   
function nf_elemDelete(ptr::Ptr{number}, cf::coeffs)
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function nf_elemCopy(a::number, cf::coeffs)
   n = julia(a)::Nemo.nf_elem
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function nf_elemGreaterZero(a::number, cf::coeffs)
   return Cint(1)
end

function nf_elemCoeffWrite(cf::coeffs, d::Cint)
   r = julia(cf)::Nemo.AnticNumberField
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function nf_elemWrite(a::number, cf::coeffs)
   n = julia(a)::Nemo.nf_elem
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

function nf_elemNeg(a::number, cf::coeffs)
   n = julia(a)::Nemo.nf_elem
   return number(-n)
end

function nf_elemInpNeg(a::number, cf::coeffs)
   R = julia(cf)::Nemo.AnticNumberField
   n = julia(a)::Nemo.nf_elem
   sub!(n, 0, n)
   return number(n, false)
end

function nf_elemInvers(a::number, cf::coeffs)
   n = julia(a)::Nemo.nf_elem
   return number(inv(n))
end

function nf_elemMult(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(n1*n2)
end

function nf_elemInpMult(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.nf_elem
   bb = julia(b)::Nemo.nf_elem
   cc = mul!(aa, aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function nf_elemAdd(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(n1 + n2)
end

function nf_elemInpAdd(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.nf_elem
   bb = julia(b)::Nemo.nf_elem
   cc = addeq!(aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function nf_elemSub(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(n1 - n2)
end

function nf_elemDiv(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function nf_elemGreater(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return Cint(n1 != n2)
end

function nf_elemEqual(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return Cint(n1 == n2)
end

function nf_elemIsZero(a::number, cf::coeffs)
   n = julia(a)::Nemo.nf_elem
   return Cint(iszero(n))
end

function nf_elemIsOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.nf_elem
   return Cint(isone(n))
end

function nf_elemIsMOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.nf_elem
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function nf_elemGcd(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.nf_elem
   n2 = julia(b)::Nemo.nf_elem
   return number(gcd(n1, n2))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function nf_elemInt(ptr::Ptr{number}, cf::coeffs)
   return Clong(0)
end

function nf_elemMPZ(b::BigInt, ptr::Ptr{number}, cf::coeffs)
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function nf_elemInitChar(cf::coeffs, p::Ptr{Void})
        
    pInit = cfunction(nf_elemInit, number, (Clong, coeffs))
    pInt = cfunction(nf_elemInt, Clong, (Ptr{number}, coeffs))
    pMPZ = cfunction(nf_elemMPZ, Void, (BigInt, Ptr{number}, coeffs))
    pInpNeg = cfunction(nf_elemInpNeg, number, (number, coeffs))
    pCopy = cfunction(nf_elemCopy, number, (number, coeffs))
    pDelete = cfunction(nf_elemDelete, Void, (Ptr{number}, coeffs))
    pAdd = cfunction(nf_elemAdd, number, (number, number, coeffs))
    pInpAdd = cfunction(nf_elemInpAdd, Void, (Ptr{number}, number, coeffs))
    pSub = cfunction(nf_elemSub, number, (number, number, coeffs))
    pMult = cfunction(nf_elemMult, number, (number, number, coeffs))
    pInpMult = cfunction(nf_elemInpMult, Void, (Ptr{number}, number, coeffs))
    pDiv = cfunction(nf_elemDiv, number, (number, number, coeffs))
    pInvers = cfunction(nf_elemInvers, number, (number, coeffs))
    pGcd = cfunction(nf_elemGcd, number, (number, number, coeffs))
    pGreater = cfunction(nf_elemGreater, Cint, (number, number, coeffs))
    pEqual = cfunction(nf_elemEqual, Cint, (number, number, coeffs))
    pIsZero = cfunction(nf_elemIsZero, Cint, (number, coeffs))
    pIsOne = cfunction(nf_elemIsOne, Cint, (number, coeffs))
    pIsMOne = cfunction(nf_elemIsMOne, Cint, (number, coeffs))
    pGreaterZero = cfunction(nf_elemGreaterZero, Cint, (number, coeffs))
    pWrite = cfunction(nf_elemWrite, Void, (number, coeffs))
    pCoeffWrite = cfunction(nf_elemCoeffWrite, Void, (coeffs, Cint))

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

function register(R::AnticNumberField)
   c = cfunction(nf_elemInitChar, Cint, (coeffs, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end

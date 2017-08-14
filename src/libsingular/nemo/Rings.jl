###############################################################################
#
#   Memory management
#
###############################################################################

function nemoRingInit(i::Clong, cf::coeffs)
   R = julia(cf)
   return number(R(i))
end
   
function nemoRingDelete(ptr::Ptr{number}, cf::coeffs)
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function nemoRingCopy(a::number, cf::coeffs)
   n = julia(a)
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function nemoRingGreaterZero(a::number, cf::coeffs)
   return Cint(1)
end

function nemoRingCoeffWrite(cf::coeffs, d::Cint)
   r = julia(cf)
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function nemoRingWrite(a::number, cf::coeffs)
   n = julia(a)
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

function nemoRingNeg(a::number, cf::coeffs)
   n = julia(a)
   return number(-n)
end

function nemoRingInpNeg(a::number, cf::coeffs)
   R = julia(cf)
   n = julia(a)
   mone = R(-1)
   mul!(n, n, mone)
   return number(n, false)
end

function nemoRingInvers(a::number, cf::coeffs)
   n = julia(a)
   return number(divexact(1, n))
end

function nemoRingMult(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(n1*n2)
end

function nemoRingInpMult(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   aa = mul!(aa, aa, bb)
   n = number(aa, false)
   nothing
end

function nemoRingAdd(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 + n2)
end

function nemoRingInpAdd(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   aa = addeq!(aa, bb)
   n = number(aa, false)
   nothing
end

function nemoRingSub(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 - n2)
end

function nemoRingDiv(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function nemoRingGreater(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 != n2)
end

function nemoRingEqual(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 == n2)
end

function nemoRingIsZero(a::number, cf::coeffs)
   n = julia(a)
   return Cint(iszero(n))
end

function nemoRingIsOne(a::number, cf::coeffs)
   n = julia(a)
   return Cint(isone(n))
end

function nemoRingIsMOne(a::number, cf::coeffs)
   n = julia(a)
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function nemoRingGcd(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(gcd(n1, n2))
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function nemoRingExtGcd(a::number, b::number, s::Ptr{number}, t::Ptr{number}, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   s1 = unsafe_load(s)
   if s1 != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(s1))
   end
   t1 = unsafe_load(t)
   if t1 != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(t1))
   end
   g1, s1, t1 = gcdx(n1, n2)
   libSingular.setindex!(s, number(s1))
   libSingular.setindex!(t, number(t1))
   return number(g1)
end

###############################################################################
#
#   Conversion
#
###############################################################################

function nemoRingInt(ptr::Ptr{number}, cf::coeffs)
   return Clong(0)
end

function nemoRingMPZ(b::BigInt, ptr::Ptr{number}, cf::coeffs)
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function nemoRingInitChar(cf::coeffs, p::Ptr{Void})
        
    pInit = cfunction(nemoRingInit, number, (Clong, coeffs))
    pInt = cfunction(nemoRingInt, Clong, (Ptr{number}, coeffs))
    pMPZ = cfunction(nemoRingMPZ, Void, (BigInt, Ptr{number}, coeffs))
    pInpNeg = cfunction(nemoRingInpNeg, number, (number, coeffs))
    pCopy = cfunction(nemoRingCopy, number, (number, coeffs))
    pDelete = cfunction(nemoRingDelete, Void, (Ptr{number}, coeffs))
    pAdd = cfunction(nemoRingAdd, number, (number, number, coeffs))
    pInpAdd = cfunction(nemoRingInpAdd, Void, (Ptr{number}, number, coeffs))
    pSub = cfunction(nemoRingSub, number, (number, number, coeffs))
    pMult = cfunction(nemoRingMult, number, (number, number, coeffs))
    pInpMult = cfunction(nemoRingInpMult, Void, (Ptr{number}, number, coeffs))
    pDiv = cfunction(nemoRingDiv, number, (number, number, coeffs))
    pInvers = cfunction(nemoRingInvers, number, (number, coeffs))
    pGcd = cfunction(nemoRingGcd, number, (number, number, coeffs))
    pExtGcd = cfunction(nemoFieldExtGcd, number, (number, number, Ptr{number}, Ptr{number}, coeffs))
    pGreater = cfunction(nemoRingGreater, Cint, (number, number, coeffs))
    pEqual = cfunction(nemoRingEqual, Cint, (number, number, coeffs))
    pIsZero = cfunction(nemoRingIsZero, Cint, (number, coeffs))
    pIsOne = cfunction(nemoRingIsOne, Cint, (number, coeffs))
    pIsMOne = cfunction(nemoRingIsMOne, Cint, (number, coeffs))
    pGreaterZero = cfunction(nemoRingGreaterZero, Cint, (number, coeffs))
    pWrite = cfunction(nemoRingWrite, Void, (number, coeffs))
    pCoeffWrite = cfunction(nemoRingCoeffWrite, Void, (coeffs, Cint))

    icxx""" 
      coeffs cf = (coeffs)($cf);
      cf->has_simple_Alloc = FALSE;  
      cf->has_simple_Inverse= FALSE;          
      cf->is_field  = FALSE;
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
      cf->cfExtGcd = (number (*)(number, number, number *, number *, const coeffs)) $pExtGcd;
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

function register(R::Nemo.Ring)
   c = cfunction(nemoRingInitChar, Cint, (coeffs, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end
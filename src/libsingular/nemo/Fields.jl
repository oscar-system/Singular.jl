###############################################################################
#
#   Memory management
#
###############################################################################

function nemoFieldInit(i::Clong, cf::coeffs)
   R = julia(cf)
   return number(R(i))
end
   
function nemoFieldDelete(ptr::Ptr{number}, cf::coeffs)
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function nemoFieldCopy(a::number, cf::coeffs)
   n = julia(a)
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function nemoFieldGreaterZero(a::number, cf::coeffs)
   return Cint(1)
end

function nemoFieldCoeffWrite(cf::coeffs, d::Cint)
   r = julia(cf)
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function nemoFieldWrite(a::number, cf::coeffs)
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

function nemoFieldNeg(a::number, cf::coeffs)
   n = julia(a)
   return number(-n)
end

function nemoFieldInpNeg(a::number, cf::coeffs)
   R = julia(cf)
   n = julia(a)
   mone = R(-1)
   mul!(n, n, mone)
   return number(n, false)
end

function nemoFieldInvers(a::number, cf::coeffs)
   n = julia(a)
   return number(inv(n))
end

function nemoFieldMult(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(n1*n2)
end

function nemoFieldInpMult(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   aa = mul!(aa, aa, bb)
   n = number(aa, false)
   nothing
end

function nemoFieldAdd(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 + n2)
end

function nemoFieldInpAdd(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   aa = addeq!(aa, bb)
   n = number(aa, false)
   nothing
end

function nemoFieldSub(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 - n2)
end

function nemoFieldDiv(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return number(divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function nemoFieldGreater(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 != n2)
end

function nemoFieldEqual(a::number, b::number, cf::coeffs)
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 == n2)
end

function nemoFieldIsZero(a::number, cf::coeffs)
   n = julia(a)
   return Cint(iszero(n))
end

function nemoFieldIsOne(a::number, cf::coeffs)
   n = julia(a)
   return Cint(isone(n))
end

function nemoFieldIsMOne(a::number, cf::coeffs)
   n = julia(a)
   return Cint(n == -1)
end

###############################################################################
#
#   Conversion
#
###############################################################################

function nemoFieldInt(ptr::Ptr{number}, cf::coeffs)
   return Clong(0)
end

function nemoFieldMPZ(b::BigInt, ptr::Ptr{number}, cf::coeffs)
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function nemoFieldInitChar(cf::coeffs, p::Ptr{Void})
        
    pInit = cfunction(nemoFieldInit, number, (Clong, coeffs))
    pInt = cfunction(nemoFieldInt, Clong, (Ptr{number}, coeffs))
    pMPZ = cfunction(nemoFieldMPZ, Void, (BigInt, Ptr{number}, coeffs))
    pInpNeg = cfunction(nemoFieldInpNeg, number, (number, coeffs))
    pCopy = cfunction(nemoFieldCopy, number, (number, coeffs))
    pDelete = cfunction(nemoFieldDelete, Void, (Ptr{number}, coeffs))
    pAdd = cfunction(nemoFieldAdd, number, (number, number, coeffs))
    pInpAdd = cfunction(nemoFieldInpAdd, Void, (Ptr{number}, number, coeffs))
    pSub = cfunction(nemoFieldSub, number, (number, number, coeffs))
    pMult = cfunction(nemoFieldMult, number, (number, number, coeffs))
    pInpMult = cfunction(nemoFieldInpMult, Void, (Ptr{number}, number, coeffs))
    pDiv = cfunction(nemoFieldDiv, number, (number, number, coeffs))
    pInvers = cfunction(nemoFieldInvers, number, (number, coeffs))
    pGreater = cfunction(nemoFieldGreater, Cint, (number, number, coeffs))
    pEqual = cfunction(nemoFieldEqual, Cint, (number, number, coeffs))
    pIsZero = cfunction(nemoFieldIsZero, Cint, (number, coeffs))
    pIsOne = cfunction(nemoFieldIsOne, Cint, (number, coeffs))
    pIsMOne = cfunction(nemoFieldIsMOne, Cint, (number, coeffs))
    pGreaterZero = cfunction(nemoFieldGreaterZero, Cint, (number, coeffs))
    pWrite = cfunction(nemoFieldWrite, Void, (number, coeffs))
    pCoeffWrite = cfunction(nemoFieldCoeffWrite, Void, (coeffs, Cint))

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

function register(R::Nemo.Field)
   c = cfunction(nemoFieldInitChar, Cint, (coeffs, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end
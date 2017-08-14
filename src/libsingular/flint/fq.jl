###############################################################################
#
#   Memory management
#
###############################################################################

function fqInit(i::Clong, cf::coeffs)
   R = julia(cf)::Nemo.FqFiniteField
   return number(R(Int(i)))
end
   
function fqDelete(ptr::Ptr{number}, cf::coeffs)
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function fqCopy(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fqGreaterZero(a::number, cf::coeffs)
   return Cint(1)
end

function fqCoeffWrite(cf::coeffs, d::Cint)
   r = julia(cf)::Nemo.FqFiniteField
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function fqWrite(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq
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

function fqNeg(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq
   return number(-n)
end

function fqInpNeg(a::number, cf::coeffs)
   R = julia(cf)::Nemo.FqFiniteField
   n = julia(a)::Nemo.fq
   ccall((:fq_neg, :libflint), Void, (Ptr{Nemo.fq}, Ptr{Nemo.fq}, Ptr{FqFiniteField}), &n, &n, &R)
   return number(n, false)
end

function fqInvers(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq
   return number(inv(n))
end

function fqMult(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1*n2)
end

function fqInpMult(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq
   bb = julia(b)::Nemo.fq
   aa = mul!(aa, aa, bb)
   n = number(aa, false)
   nothing
end

function fqAdd(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1 + n2)
end

function fqInpAdd(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq
   bb = julia(b)::Nemo.fq
   aa = addeq!(aa, bb)
   n = number(aa, false)
   nothing
end

function fqSub(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1 - n2)
end

function fqDiv(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fqGreater(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return Cint(n1 != n2)
end

function fqEqual(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return Cint(n1 == n2)
end

function fqIsZero(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq
   return Cint(iszero(n))
end

function fqIsOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq
   return Cint(isone(n))
end

function fqIsMOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.fq
   return Cint(n == -1)
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fqInt(ptr::Ptr{number}, cf::coeffs)
   return Clong(0)
end

function fqMPZ(b::BigInt, ptr::Ptr{number}, cf::coeffs)
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function fqInitChar(cf::coeffs, p::Ptr{Void})
        
    pInit = cfunction(fqInit, number, (Clong, coeffs))
    pInt = cfunction(fqInt, Clong, (Ptr{number}, coeffs))
    pMPZ = cfunction(fqMPZ, Void, (BigInt, Ptr{number}, coeffs))
    pInpNeg = cfunction(fqInpNeg, number, (number, coeffs))
    pCopy = cfunction(fqCopy, number, (number, coeffs))
    pDelete = cfunction(fqDelete, Void, (Ptr{number}, coeffs))
    pAdd = cfunction(fqAdd, number, (number, number, coeffs))
    pInpAdd = cfunction(fqInpAdd, Void, (Ptr{number}, number, coeffs))
    pSub = cfunction(fqSub, number, (number, number, coeffs))
    pMult = cfunction(fqMult, number, (number, number, coeffs))
    pInpMult = cfunction(fqInpMult, Void, (Ptr{number}, number, coeffs))
    pDiv = cfunction(fqDiv, number, (number, number, coeffs))
    pInvers = cfunction(fqInvers, number, (number, coeffs))
    pGreater = cfunction(fqGreater, Cint, (number, number, coeffs))
    pEqual = cfunction(fqEqual, Cint, (number, number, coeffs))
    pIsZero = cfunction(fqIsZero, Cint, (number, coeffs))
    pIsOne = cfunction(fqIsOne, Cint, (number, coeffs))
    pIsMOne = cfunction(fqIsMOne, Cint, (number, coeffs))
    pGreaterZero = cfunction(fqGreaterZero, Cint, (number, coeffs))
    pWrite = cfunction(fqWrite, Void, (number, coeffs))
    pCoeffWrite = cfunction(fqCoeffWrite, Void, (coeffs, Cint))

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

function register(R::FqFiniteField)
   c = cfunction(fqInitChar, Cint, (coeffs, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end
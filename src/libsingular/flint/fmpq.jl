###############################################################################
#
#   Memory management
#
###############################################################################

function fmpqInit(i::Clong, cf::coeffs)
   return number(Nemo.fmpq(i))
end
   
function fmpqDelete(ptr::Ptr{number}, cf::coeffs)
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function fmpqCopy(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fmpqGreaterZero(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   return Cint(n > 0)
end

function fmpqCoeffWrite(cf::coeffs, d::Cint)
   r = julia(cf)::Nemo.FlintRationalField
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function fmpqWrite(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
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

function fmpqNeg(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   return number(-n)
end

function fmpqInpNeg(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   ccall((:fmpq_neg, :libflint), Void, (Ptr{Nemo.fmpq}, Ptr{Nemo.fmpq}), &n, &n)
   return number(n, false)
end

function fmpqInvers(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   return number(inv(n))
end

function fmpqMult(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fmpq
   n2 = julia(b)::Nemo.fmpq
   return number(n1*n2)
end

function fmpqInpMult(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fmpq
   bb = julia(b)::Nemo.fmpq
   cc = mul!(aa, aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function fmpqAdd(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fmpq
   n2 = julia(b)::Nemo.fmpq
   return number(n1 + n2)
end

function fmpqInpAdd(a::Ptr{number}, b::number, cf::coeffs)
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fmpq
   bb = julia(b)::Nemo.fmpq
   cc = addeq!(aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function fmpqSub(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fmpq
   n2 = julia(b)::Nemo.fmpq
   return number(n1 - n2)
end

function fmpqDiv(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fmpq
   n2 = julia(b)::Nemo.fmpq
   return number(divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fmpqGreater(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fmpq
   n2 = julia(b)::Nemo.fmpq
   return Cint(n1 > n2)
end

function fmpqEqual(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fmpq
   n2 = julia(b)::Nemo.fmpq
   return Cint(n1 == n2)
end

function fmpqIsZero(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   return Cint(iszero(n))
end

function fmpqIsOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   return Cint(isone(n))
end

function fmpqIsMOne(a::number, cf::coeffs)
   n = julia(a)::Nemo.fmpq
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function fmpqGcd(a::number, b::number, cf::coeffs)
   n1 = julia(a)::Nemo.fmpq
   n2 = julia(b)::Nemo.fmpq
   return number(gcd(n1, n2))
end

###############################################################################
#
#   Subring GCD
#
###############################################################################

function fmpqSubringGcd(a::number, b::number, cf::coeffs)
   n1 = numerator(julia(a)::fmpq)
   n2 = numerator(julia(b)::fmpq)
   return number(fmpq(gcd(n1, n2)))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fmpqInt(ptr::Ptr{number}, cf::coeffs)
   n = julia(unsafe_load(ptr))::fmpq
   return Clong(num)
end

function fmpqMPZ(b::BigInt, ptr::Ptr{number}, cf::coeffs)
   n = julia(unsafe_load(ptr))::fmpq
   z = convert(BigInt, num(n))
   bptr = pointer_from_objref(b)
   zptr = pointer_from_objref(z)
   icxx"""mpz_init_set((__mpz_struct *) $bptr, (mpz_ptr) $zptr);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function fmpqInitChar(cf::coeffs, p::Ptr{Void})
        
    pInit = cfunction(fmpqInit, number, (Clong, coeffs))
    pInt = cfunction(fmpqInt, Clong, (Ptr{number}, coeffs))
    pMPZ = cfunction(fmpqMPZ, Void, (BigInt, Ptr{number}, coeffs))
    pInpNeg = cfunction(fmpqInpNeg, number, (number, coeffs))
    pCopy = cfunction(fmpqCopy, number, (number, coeffs))
    pDelete = cfunction(fmpqDelete, Void, (Ptr{number}, coeffs))
    pAdd = cfunction(fmpqAdd, number, (number, number, coeffs))
    pInpAdd = cfunction(fmpqInpAdd, Void, (Ptr{number}, number, coeffs))
    pSub = cfunction(fmpqSub, number, (number, number, coeffs))
    pMult = cfunction(fmpqMult, number, (number, number, coeffs))
    pInpMult = cfunction(fmpqInpMult, Void, (Ptr{number}, number, coeffs))
    pDiv = cfunction(fmpqDiv, number, (number, number, coeffs))
    pInvers = cfunction(fmpqInvers, number, (number, coeffs))
    pGcd = cfunction(fmpqGcd, number, (number, number, coeffs))
    pSubringGcd = cfunction(fmpqSubringGcd, number, (number, number, coeffs))
    pGreater = cfunction(fmpqGreater, Cint, (number, number, coeffs))
    pEqual = cfunction(fmpqEqual, Cint, (number, number, coeffs))
    pIsZero = cfunction(fmpqIsZero, Cint, (number, coeffs))
    pIsOne = cfunction(fmpqIsOne, Cint, (number, coeffs))
    pIsMOne = cfunction(fmpqIsMOne, Cint, (number, coeffs))
    pGreaterZero = cfunction(fmpqGreaterZero, Cint, (number, coeffs))
    pWrite = cfunction(fmpqWrite, Void, (number, coeffs))
    pCoeffWrite = cfunction(fmpqCoeffWrite, Void, (coeffs, Cint))

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
      cf->cfSubringGcd = (numberfunc) $pSubringGcd;
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

function register(R::FlintRationalField)
   c = cfunction(fmpqInitChar, Cint, (coeffs, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end

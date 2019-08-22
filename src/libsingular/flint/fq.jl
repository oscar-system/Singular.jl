###############################################################################
#
#   Memory management
#
###############################################################################

function fqInit(i::Clong, cf::Ptr{Cvoid})
   R = julia(cf)::Nemo.FqFiniteField
   return number(R(Int(i)))
end
   
function fqDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function fqCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function fqGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   return Cint(1)
end

function fqCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
   r = julia(cf)::Nemo.FqFiniteField
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function fqWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
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

function fqNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return number(-n)
end

function fqInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   R = julia(cf)::Nemo.FqFiniteField
   n = julia(a)::Nemo.fq
   ccall((:fq_neg, :libflint), Void, (Ptr{Nemo.fq}, Ptr{Nemo.fq}, Ptr{FqFiniteField}), &n, &n, &R)
   return number(n, false)
end

function fqInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return number(inv(n))
end

function fqMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1*n2)
end

function fqInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq
   bb = julia(b)::Nemo.fq
   cc = mul!(aa, aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function fqAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1 + n2)
end

function fqInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)::Nemo.fq
   bb = julia(b)::Nemo.fq
   cc = addeq!(aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function fqSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(n1 - n2)
end

function fqDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(divexact(n1, n2))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function fqGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return Cint(n1 != n2)
end

function fqEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return Cint(n1 == n2)
end

function fqIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return Cint(iszero(n))
end

function fqIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return Cint(isone(n))
end

function fqIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)::Nemo.fq
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function fqGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)::Nemo.fq
   n2 = julia(b)::Nemo.fq
   return number(gcd(n1, n2))
end

###############################################################################
#
#   Conversion
#
###############################################################################

function fqInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   return Clong(0)
end

function fqMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function fqInitChar(cf::Ptr{Cvoid}, p::Ptr{Void})
        
    pInit = cfunction(fqInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    pInt = cfunction(fqInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pMPZ = cfunction(fqMPZ, Void, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pInpNeg = cfunction(fqInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pCopy = cfunction(fqCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pDelete = cfunction(fqDelete, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pAdd = cfunction(fqAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInpAdd = cfunction(fqInpAdd, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    pSub = cfunction(fqSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pMult = cfunction(fqMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInpMult = cfunction(fqInpMult, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    pDiv = cfunction(fqDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInvers = cfunction(fqInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pGcd = cfunction(fqGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pGreater = cfunction(fqGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pEqual = cfunction(fqEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pIsZero = cfunction(fqIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pIsOne = cfunction(fqIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pIsMOne = cfunction(fqIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pGreaterZero = cfunction(fqGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pWrite = cfunction(fqWrite, Void, (Ptr{Cvoid}, Ptr{Cvoid}))
    pCoeffWrite = cfunction(fqCoeffWrite, Void, (Ptr{Cvoid}, Cint))

    return Cint(0)
end

function register(R::FqFiniteField)
   c = cfunction(fqInitChar, Cint, (Ptr{Cvoid}, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end

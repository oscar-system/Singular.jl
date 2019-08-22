###############################################################################
#
#   Memory management
#
###############################################################################

function nemoRingInit(i::Clong, cf::Ptr{Cvoid})
   R = julia(cf)
   return number(R(i))
end
   
function nemoRingDelete(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n = unsafe_load(ptr)
   if n != C_NULL
      number_pop!(nemoNumberID, Ptr{Void}(n))
   end
   nothing
end

function nemoRingCopy(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(deepcopy(n))
end

###############################################################################
#
#   Printing
#
###############################################################################

function nemoRingGreaterZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   return Cint(1)
end

function nemoRingCoeffWrite(cf::Ptr{Cvoid}, d::Cint)
   r = julia(cf)
   str = string(r)
   icxx"""PrintS($str);"""
   nothing
end

function nemoRingWrite(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
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

function nemoRingNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(-n)
end

function nemoRingInpNeg(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   R = julia(cf)
   n = julia(a)
   mone = R(-1)
   n_new = mul!(n, n, mone)
   return number(n_new, n)
end

function nemoRingInvers(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return number(divexact(1, n))
end

function nemoRingMult(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1*n2)
end

function nemoRingInpMult(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   cc = mul!(aa, aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function nemoRingAdd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 + n2)
end

function nemoRingInpAdd(a::Ptr{Ptr{Cvoid}}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   r = unsafe_load(a)
   aa = julia(r)
   bb = julia(b)
   cc = addeq!(aa, bb)
   n = number(cc, aa)
   unsafe_store!(a, n, 1)
   nothing
end

function nemoRingSub(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(n1 - n2)
end

function nemoRingDiv(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(divexact(n1, n2))
end

function nemoRingDivBy(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return Cint(divides(n1, n2)[1])
end

###############################################################################
#
#   Comparison
#
###############################################################################

function nemoRingGreater(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 != n2)
end

function nemoRingEqual(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return Cint(n1 == n2)
end

function nemoRingIsZero(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(iszero(n))
end

function nemoRingIsOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(isone(n))
end

function nemoRingIsMOne(a::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n = julia(a)
   return Cint(n == -1)
end

###############################################################################
#
#   GCD
#
###############################################################################

function nemoRingGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
   return number(gcd(n1, n2))
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function nemoRingExtGcd(a::Ptr{Cvoid}, b::Ptr{Cvoid}, s::Ptr{Ptr{Cvoid}}, t::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   n1 = julia(a)
   n2 = julia(b)
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

function nemoRingInt(ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   return Clong(0)
end

function nemoRingMPZ(b::BigInt, ptr::Ptr{Ptr{Cvoid}}, cf::Ptr{Cvoid})
   bptr = pointer_from_objref(b)
   icxx"""mpz_init_set_si((__mpz_struct *) $bptr, 0);"""
   nothing
end

###############################################################################
#
#   InitChar
#
###############################################################################

function nemoRingInitChar(cf::Ptr{Cvoid}, p::Ptr{Void})
        
    pInit = cfunction(nemoRingInit, Ptr{Cvoid}, (Clong, Ptr{Cvoid}))
    pInt = cfunction(nemoRingInt, Clong, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pMPZ = cfunction(nemoRingMPZ, Void, (BigInt, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pInpNeg = cfunction(nemoRingInpNeg, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pCopy = cfunction(nemoRingCopy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pDelete = cfunction(nemoRingDelete, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pAdd = cfunction(nemoRingAdd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInpAdd = cfunction(nemoRingInpAdd, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    pSub = cfunction(nemoRingSub, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pMult = cfunction(nemoRingMult, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInpMult = cfunction(nemoRingInpMult, Void, (Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Cvoid}))
    pDiv = cfunction(nemoRingDiv, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pDivBy = cfunction(nemoRingDivBy, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pInvers = cfunction(nemoRingInvers, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}))
    pGcd = cfunction(nemoRingGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pExtGcd = cfunction(nemoRingExtGcd, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}))
    pGreater = cfunction(nemoRingGreater, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pEqual = cfunction(nemoRingEqual, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
    pIsZero = cfunction(nemoRingIsZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pIsOne = cfunction(nemoRingIsOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pIsMOne = cfunction(nemoRingIsMOne, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pGreaterZero = cfunction(nemoRingGreaterZero, Cint, (Ptr{Cvoid}, Ptr{Cvoid}))
    pWrite = cfunction(nemoRingWrite, Void, (Ptr{Cvoid}, Ptr{Cvoid}))
    pCoeffWrite = cfunction(nemoRingCoeffWrite, Void, (Ptr{Cvoid}, Cint))

    return Cint(0)
end

function register(R::Nemo.Ring)
   c = cfunction(nemoRingInitChar, Cint, (Ptr{Cvoid}, Ptr{Void}))
   ptr = @cxx n_unknown
   return nRegister(ptr, c)
end

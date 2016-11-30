# initialise a coefficient ring
function nInitChar(n::n_coeffType, p::Ptr{Void})
   return icxx"""nInitChar($n, $p);"""
end

# kill a coefficient ring
function nKillChar(cf::coeffs)
   icxx"""nKillChar($cf);"""
end

function nCoeff_has_simple_Alloc(cf::coeffs)
   icxx"""nCoeff_has_simple_Alloc($cf);""" > 0
end

# initialise a number from an Int
function n_Init(i::Int, cf::coeffs) 
   icxx"""n_Init($i, $cf);"""
end

# initialise a number from an Int
function n_Copy(n::number, cf::coeffs) 
   icxx"""n_Copy($n, $cf);"""
end

# initialise a number from an mpz
function n_InitMPZ(b::BigInt, cf::coeffs)
    bb = __mpz_struct(pointer_from_objref(b))
    icxx"""n_InitMPZ($bb, $cf);"""
end

# delete a number
function n_Delete(n::number, cf::coeffs) 
   icxx"""number t = $n; if (t != NULL) n_Delete(&t, $cf);"""
end

# write a number to a Singular string
function n_Write(n::number_ref, cf::coeffs, bShortOut::Bool = false)
   d = Int(bShortOut)
   icxx"""n_Write($n, $cf, $d);"""
end

function n_Add(a::number, b::number, cf::coeffs)
   icxx"""n_Add($a, $b, $cf);"""
end

function n_Sub(a::number, b::number, cf::coeffs)
   icxx"""n_Sub($a, $b, $cf);"""
end

function n_Mult(a::number, b::number, cf::coeffs)
   icxx"""n_Mult($a, $b, $cf);"""
end

function n_Neg(n::number, cf::coeffs)
   icxx"""number nn = n_Copy($n, $cf); nn = n_InpNeg(nn, $cf); nn;"""
end

function n_Invers(a::number, cf::coeffs)
   icxx"""n_Invers($a, $cf);"""
end

function n_ExactDiv(a::number, b::number, cf::coeffs)
   icxx"""n_ExactDiv($a, $b, $cf);"""
end

function n_Div(a::number, b::number, cf::coeffs)
   icxx"""number z = n_Div($a, $b, $cf); n_Normalize(z, $cf); z;"""
end

function n_Power(a::number, b::Int, cf::coeffs)
   icxx"""number res; n_Power($a, $b, &res, $cf); res;"""
end

function n_Gcd(a::number, b::number, cf::coeffs)
   icxx"""n_Gcd($a, $b, $cf);"""
end

function n_Lcm(a::number, b::number, cf::coeffs)
   icxx"""n_Lcm($a, $b, $cf);"""
end

function n_IsZero(a::number, cf::coeffs)
   icxx"""n_IsZero($a, $cf);""" > 0
end

function n_IsOne(a::number, cf::coeffs)
   icxx"""n_IsOne($a, $cf);""" > 0
end

function n_Greater(a::number, b::number, cf::coeffs)
   icxx"""n_Greater($a, $b, $cf);""" > 0
end

function n_Equal(a::number, b::number, cf::coeffs)
   icxx"""n_Equal($a, $b, $cf);""" > 0
end

function n_InpAdd(a::number_ref, b::number, cf::coeffs)
   icxx"""n_InpAdd($a, $b, $cf);"""
end

function n_InpMult(a::number_ref, b::number, cf::coeffs)
   icxx"""n_InpMult($a, $b, $cf);"""
end

function n_QuotRem(a::number, b::number, cf::coeffs)
   q = n_Init(0, cf)
   r = n_Init(0, cf)
   icxx"""number qq = n_Init(0, $cf); $r = n_QuotRem($a, $b, &qq, $cf); $q = qq;"""
   return q[], r[]
end

function n_Rem(a::number, b::number, cf::coeffs)
   q = n_Init(0, cf)
   r = n_Init(0, cf)
   icxx"""number qq = n_Init(0, $cf); $r = n_QuotRem($a, $b, &qq, $cf); $q = qq;"""
   return r[]
end

# create a Singular string environment
function StringSetS(m) 
   @cxx StringSetS(pointer(m))
end

# end a Singular string environment
function StringEndS() 
   icxx"""StringEndS();"""
end

# omalloc free
function omFree{T}(m :: Ptr{T})
   icxx"""omFree((void*) $m);"""
end
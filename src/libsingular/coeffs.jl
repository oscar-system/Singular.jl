# return a function to convert between rings
function n_SetMap(src::coeffs, dst::coeffs)
   icxx"""n_SetMap($src, $dst);"""
end

function nApplyMapFunc(f::Cxx.CppFptr, n::number, src::coeffs, dst::coeffs)
   icxx"""$f($n, $src, $dst);"""
end

function nCoeff_has_simple_Alloc(cf::coeffs)
   icxx"""nCoeff_has_simple_Alloc($cf);""" > 0
end

# initialise a number from an mpz
function n_InitMPZ(b::BigInt, cf::coeffs)
    bb = __mpz_struct(pointer_from_objref(b))
    return n_InitMPZ_internal(bb,cf)
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

function n_GetNumerator(a::number_ref, cf::coeffs)
   icxx"""return n_GetNumerator($a, $cf);"""
end

function n_GetDenom(a::number_ref, cf::coeffs)
   icxx"""return n_GetDenom($a, $cf);"""
end

function n_Power(a::number, b::Int, cf::coeffs)
   icxx"""number res; n_Power($a, $b, &res, $cf); res;"""
end

function n_Gcd(a::number, b::number, cf::coeffs)
   icxx"""n_Gcd($a, $b, $cf);"""
end

function n_SubringGcd(a::number, b::number, cf::coeffs)
   icxx"""n_SubringGcd($a, $b, $cf);"""
end

function n_Lcm(a::number, b::number, cf::coeffs)
   icxx"""n_Lcm($a, $b, $cf);"""
end

function n_ExtGcd(a::number, b::number, s::Ptr{number}, t::Ptr{number}, cf:: coeffs)
   icxx"""n_ExtGcd($a, $b, $s, $t, $cf);"""
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

function n_GreaterZero(a::number, cf::coeffs)
   icxx"""n_GreaterZero($a, $cf);""" > 0
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

function n_QuotRem(a::number, b::number, p::Ptr{number}, cf::coeffs)
   icxx"""n_QuotRem($a, $b, $p, $cf);"""
end

function n_Rem(a::number, b::number, cf::coeffs)
   icxx"""number qq; n_QuotRem($a, $b, &qq, $cf);"""
end

function n_IntMod(a::number, b::number, cf::coeffs)
   icxx"""n_IntMod($a, $b, $cf);"""
end

function n_Farey(a::number, b::number, cf::coeffs)
   icxx"""n_Farey($a, $b, $cf);"""
end

function n_ChineseRemainderSym(a::Array{number, 1}, b::Array{number, 1}, n::Cint, signed::Cint, cf::coeffs)
   p1 = pointer(a)
   p2 = pointer(b)
   icxx"""CFArray inv_cache($n); n_ChineseRemainderSym($p1, $p2, $n, $signed, inv_cache, $cf);"""
end

function n_Param(a::Cint, cf::coeffs)
   icxx"""n_Param($a, $cf);"""
end

# create a Singular string environment
function StringSetS(m) 
   @cxx StringSetS(pointer(m))
end

# end a Singular string environment
function StringEndS() 
   icxx"""StringEndS();"""
end

function omAlloc0(size :: Csize_t)
   icxx"""(void*) omAlloc0($size);"""
end

# omalloc free
function omFree{T}(m :: Ptr{T})
   icxx"""omFree((void*) $m);"""
end

###############################################################################
#
#   Setting a Ptr{number} to a number
#
###############################################################################

function setindex!(ptr::Ptr{number}, n::number)
   icxx"""*$ptr = $n;"""
   nothing
end

###############################################################################
#
#   Conversion between numbers that wrap jl_value_t's and julia values
#
###############################################################################

type live_cache
   num::Int
   A::Array{Nemo.RingElem, 1}
end

function nRegister(t::n_coeffType, f::Ptr{Void})
   return icxx"""return nRegister($t, (cfInitCharProc)$f);"""
end

const nemoNumberID = Base.Dict{UInt, live_cache}()

function julia(cf::coeffs)
   ptr = @cxx cf->data
   return unsafe_pointer_to_objref(ptr)
end

function julia(p::number)
    ptr = icxx"""return (void*)$p;"""
    return unsafe_pointer_to_objref(ptr)
end

function number_pop!(D::Base.Dict{UInt, live_cache}, ptr::Ptr{Void})
   iptr = reinterpret(UInt, ptr) >> 24
   val = D[iptr]
   val.num -= 1
   if val.num == 0
      pop!(D, iptr)
   end
end

function number{T <: RingElem}(j::T, cache::Bool=true)
    ptr = pointer_from_objref(j)
    n = number(ptr)
    if cache
       iptr = reinterpret(UInt, ptr) >> 24
       if !haskey(nemoNumberID, iptr)
          nemoNumberID[iptr] = live_cache(0, Array{Nemo.RingElem}(64))
       end
       val = nemoNumberID[iptr]
       val.num += 1
       push!(val.A, j)
    end
    return n
end

function number{T <: RingElem}(j::T, j_old::T)
    if j == j_old   # for inplace operations
        return number(j, false)
    else
        number_pop!(nemoNumberID, pointer_from_objref(j_old))
        return number(j)
    end
end

include("flint/fmpz.jl")
include("flint/fmpq.jl")
include("flint/fq.jl")
include("flint/fq_nmod.jl")
include("antic/nf_elem.jl")
include("nemo/Rings.jl")
include("nemo/Fields.jl")


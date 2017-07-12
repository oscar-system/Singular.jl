function rDefault(cf::coeffs, vars::Array{Ptr{UInt8}, 1}, ord::rRingOrder_t)
   len = length(vars)
   ptr = pointer(vars)
   r = icxx"""rDefault($cf, $len, $ptr, $ord);"""
   icxx"""$r->ShortOut = 0;"""
   return r
end

function rDefault{T}(cf::coeffs, vars::Array{T,1}, ord::Array{rRingOrder_t, 1}, 
   blk0::Array{Cint, 1}, blk1::Array{Cint, 1})
   wvhdl = Ptr{Ptr{Cint}}(C_NULL)
   len = length(vars)
   ptr = Ptr{Ptr{Cuchar}}(pointer(vars))
   ordlen = length(ord)
   ordptr =  pointer(ord)
   blk0ptr = pointer(blk0)
   blk1ptr = pointer(blk1)
   r = icxx"""rDefault($cf, $len, $ptr, $ordlen, $ordptr, $blk0ptr, $blk1ptr, $wvhdl);"""
   icxx"""$r->ShortOut = 0;"""
   return r
end

function rDelete(r::ring)
   icxx"""rDelete($r);"""
end

function rString(r::ring)
   icxx"""rString($r);"""
end

function rChar(r::ring)
   icxx"""rChar($r);"""
end

function  rGetVar(i::Cint, r::ring)
   icxx"""rGetVar($i, $r);"""
end

function  rVar(r::ring)
   icxx"""rVar($r);"""
end

function p_Delete(p::poly, r::ring)
   icxx"""p_Delete(&$p, $r);"""
end

function p_Copy(p::poly, r::ring)
   icxx"""p_Copy($p, $r);"""
end

function p_IsOne(p::poly, r::ring)
   icxx"""p_IsOne($p, $r);"""
end

function p_IsUnit(p::poly, r::ring)
   icxx"""p_IsUnit($p, $r);"""
end

function p_GetExp(p::poly, i::Cint, r::ring)
   icxx"""p_GetExp($p, $i, $r);"""
end

function p_String(p::poly, r::ring)
   icxx"""p_String($p, $r);"""
end

function p_ISet(i:: Int, r::ring)
    icxx"""p_ISet($i, $r);"""
end

function p_NSet(n::number, r::ring)
   icxx"""p_NSet($n, $r);"""
end

function pLength(a::poly)
   icxx"""pLength($a);"""
end

function pNext(a::poly)
   icxx"""poly p = pNext($a); p;"""
end

function p_Neg(a::poly, r::ring)
   icxx"""p_Neg($a, $r);"""
end

function pGetCoeff(a::poly)
   icxx"""number p = pGetCoeff($a); p;"""
end

function p_Add_q(a::poly, b::poly, r::ring)
   icxx"""p_Add_q($a, $b, $r);"""
end

function p_Sub(a::poly, b::poly, r::ring)
   icxx"""p_Sub($a, $b, $r);"""
end

function p_Mult_q(a::poly, b::poly, r::ring)
   icxx"""p_Mult_q($a, $b, $r);"""
end

function p_Power(a::poly, n::Cint, r::ring)
   icxx"""p_Power($a, $n, $r);"""
end

function p_EqualPolys(a::poly, b::poly, r::ring)
   icxx"""p_EqualPolys($a, $b, $r);"""
end

function singclap_pdivide(a::poly, b::poly, r::ring)
   icxx"""singclap_pdivide($a, $b, $r);"""
end

function singclap_gcd(a::poly, b::poly, r::ring)
   icxx"""singclap_gcd($a, $b, $r);"""
end

function singclap_extgcd(a::poly, b::poly, res::poly_ref, s::poly_ref, t::poly_ref, r::ring)
   icxx"""singclap_extgcd($a, $b, $res, $s, $t, $r);"""
end

function p_Content(a::poly, r::ring)
   icxx"""p_Content($a, $r);"""
end

function p_Reduce(p::poly, G::ideal, R::ring)
   icxx"""const ring origin = currRing;
          rChangeCurrRing($R);
          poly res = kNF($G, $R->qideal, $p);
          rChangeCurrRing(origin);
          res;
       """
end

function p_Reduce(p::ideal, G::ideal, R::ring)
   icxx"""const ring origin = currRing;
          rChangeCurrRing($R);
          ideal res = kNF($G, $R->qideal, $p);
          rChangeCurrRing(origin);
          res;
       """
end

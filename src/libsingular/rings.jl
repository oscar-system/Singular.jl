function rDefault(cf::coeffs, vars::Array{Ptr{UInt8}, 1}, ord::rRingOrder_t)
   len = length(vars)
   ptr = pointer(vars)
   icxx"""rDefault($cf, $len, $ptr, $ord);"""
end

function rDelete(r::ring)
   icxx"""rDelete($r);"""
end

function p_Delete(p::poly, r::ring)
   icxx"""p_Delete(&$p, $r);"""
end

function p_Copy(p::poly, r::ring)
   icxx"""p_Copy($p, $r);"""
end

function p_ISet(i:: Int, r::ring)
    icxx"""p_ISet($i, $r);"""
end

function p_NSet(n::number, r::ring)
   icxx"""p_NSet($n, $r);"""
end

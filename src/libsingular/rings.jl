# function rDefault{T}(cf::coeffs, vars::Array{T,1}, ord::Array{rRingOrder_t, 1}, 
#    blk0::Array{Cint, 1}, blk1::Array{Cint, 1}, bitmask::Culong)
#    wvhdl = Ptr{Ptr{Cint}}(C_NULL)
#    len = length(vars)
#    ptr = Ptr{Ptr{Cuchar}}(pointer(vars))
#    ordlen = length(ord)
#    ordptr =  pointer(ord)
#    blk0ptr = pointer(blk0)
#    blk1ptr = pointer(blk1)
#    r = icxx"""rDefault($cf, $len, $ptr, $ordlen, $ordptr, $blk0ptr, $blk1ptr, $wvhdl, $bitmask);"""
#    icxx"""$r->ShortOut = 0;"""
#    return r
# end

function p_GetExpVL(p::poly, ev::Array{Clong, 1}, r::ring)
   ptr = pointer(ev)
   p_GetExpVL_internal(p, ptr, r)
end

function p_SetExpV(p::poly, ev::Array{Cint, 1}, r::ring)
   @assert ev[1] == 0
   ptr = pointer(ev)
   p_SetExpV_internal(p, ptr, r)
end

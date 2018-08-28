function rDefault{T}(cf::coeffs, vars::Array{T,1}, ord::Array{rRingOrder_t, 1}, 
   blk0::Array{Cint, 1}, blk1::Array{Cint, 1}, bitmask::Culong)
   blk0ptr = pointer(blk0)
   blk1ptr = pointer(blk1)
   r = rDefault_long_helper(cf, vars, ord, blk0ptr, blk1ptr, bitmask);
   return r
end

function p_GetExpVL(p::poly, ev::Array{Clong, 1}, r::ring)
   ptr = pointer(ev)
   p_GetExpVL_internal(p, ptr, r)
end

function p_SetExpV(p::poly, ev::Array{Cint, 1}, r::ring)
   @assert ev[1] == 0
   ptr = pointer(ev)
   p_SetExpV_internal(p, ptr, r)
end

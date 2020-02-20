function rDefault(cf::coeffs_ptr, vars::Array{T,1}, ord::Array{rRingOrder_t, 1},
   blk0::Array{Cint, 1}, blk1::Array{Cint, 1}, bitmask::Culong) where {T}
   blk0ptr = pointer(blk0)
   blk1ptr = pointer(blk1)
   r = rDefault_long_helper(cf, vars, ord, blk0ptr, blk1ptr, bitmask);
   return r
end

function p_GetExpVL(p::poly_ptr, ev::Array{Clong, 1}, r::ring_ptr)
   ptr = pointer(ev)
   p_GetExpVL_internal(p, ptr, r)
end

function p_SetExpV(p::poly_ptr, ev::Array{Cint, 1}, r::ring_ptr)
   @assert ev[1] == 0
   ptr = pointer(ev)
   p_SetExpV_internal(p, ptr, r)
end

function p_ExtGcd(a::poly_ptr, b::poly_ptr, res::Ptr{poly_ptr}, s::Ptr{poly_ptr}, t::Ptr{poly_ptr}, r::ring_ptr)
    sp = reinterpret(Ptr{Nothing},s)
    tp = reinterpret(Ptr{Nothing},t)
    rp = reinterpret(Ptr{Nothing},res)
    return p_ExtGcd_internal(a, b, rp, sp, tp, r)
 end

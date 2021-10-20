function p_GetExpVL(p::poly_ptr, ev::Vector{Clong}, r::ring_ptr)
   ptr = pointer(ev)
   p_GetExpVL_internal(p, ptr, r)
end

function p_GetExpVLV(p::poly_ptr, ev::Vector{Clong}, r::ring_ptr)
   ptr = pointer(ev)
   p_GetExpVLV_internal(p, ptr, r)
end

function p_SetExpV(p::poly_ptr, ev::Vector{Cint}, r::ring_ptr)
   @assert ev[1] == 0
   ptr = pointer(ev)
   p_SetExpV_internal(p, ptr, r)
end

function p_SetExpVL(p::poly_ptr, ev::Vector{Clong}, r::ring_ptr)
   ptr = pointer(ev)
   p_SetExpVL_internal(p, ptr, r)
end

function p_SetExpVLV(p::poly_ptr, ev::Vector{Clong}, c::Clong, r::ring_ptr)
   ptr = pointer(ev)
   p_SetExpVLV_internal(p, ptr, c, r)
end

function p_ExtGcd(a::poly_ptr, b::poly_ptr, res::Ptr{poly_ptr}, s::Ptr{poly_ptr}, t::Ptr{poly_ptr}, r::ring_ptr)
    sp = reinterpret(Ptr{Nothing},s)
    tp = reinterpret(Ptr{Nothing},t)
    rp = reinterpret(Ptr{Nothing},res)
    return p_ExtGcd_internal(a, b, rp, sp, tp, r)
end

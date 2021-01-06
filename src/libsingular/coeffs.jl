# initialise a number from an mpz
function n_InitMPZ(b::BigInt, cf::coeffs_ptr)
    bb = pointer_from_objref(b)
    return n_InitMPZ_internal(bb, cf)
end

# get an mpz from a number
function n_GetMPZ(s::number_ptr, r::coeffs_ptr)
   res = BigInt(1)
   resp = pointer_from_objref(res)
   n_GetMPZ_internal(resp, s, r)
   return res
end

# write a number to a Singular string
function n_Write(n::number_ptr, cf::coeffs_ptr, bShortOut::Bool = false)
   d = Int(bShortOut)
   n_Write_internal(n, cf, d);
end

function n_ExtGcd(a::number_ptr, b::number_ptr, s::Ptr{number_ptr}, t::Ptr{number_ptr}, cf:: coeffs_ptr)
   sp = reinterpret(Ptr{Nothing}, s)
   tp = reinterpret(Ptr{Nothing}, t)
   return n_ExtGcd_internal(a, b, sp, tp, cf);
end

function n_QuotRem(a::number_ptr, b::number_ptr, p::Ptr{number_ptr}, cf::coeffs_ptr)
   pp = reinterpret(Ptr{Nothing}, p)
   n_QuotRem_internal(a, b, pp, cf)
end

function n_ChineseRemainderSym(a::Array{number_ptr, 1}, b::Array{number_ptr, 1}, n::Cint, signed::Cint, cf::coeffs_ptr)
   p1 = reinterpret(Ptr{Nothing}, pointer(a))
   p2 = reinterpret(Ptr{Nothing}, pointer(b))
   return n_ChineseRemainderSym_internal(p1, p2, n, signed, cf)
end

# create a Singular string environment
function StringSetS(m)
   StringSetS_internal(m)
end

# omalloc free
function omFree(m :: Ptr{T}) where {T}
   mp = reinterpret(Ptr{Nothing}, m)
   omFree_internal(m)
end

###############################################################################
#
#   Setting a Ptr{number} to a number
#
###############################################################################

function setindex!(ptr::Ptr{number_ptr}, n::number_ptr)
   pp = reinterpret(Ptr{Nothing}, ptr)
   setindex_internal(pp, n)
end

###############################################################################
#
#   Conversion between numbers that wrap jl_value_t's and julia values
#
###############################################################################

import Nemo

mutable struct live_cache
   num::Int
   A::Array{Nemo.RingElem, 1}
end

const nemoNumberID = Base.Dict{UInt, live_cache}()

function julia(cf::coeffs_ptr)
   ptr = get_coeff_data(cf)
   return unsafe_pointer_to_objref(ptr)
end

function julia(p::Ptr{Cvoid})
    return unsafe_pointer_to_objref(p)
end

function number_pop!(D::Base.Dict{UInt, live_cache}, ptr::Ptr{Cvoid})
    iptr = reinterpret(UInt, ptr) >> 24
    if haskey(D, iptr)
        val = D[iptr]
        val.num -= 1
        if val.num == 0
           delete!(D, iptr)
        end
    end
end

function number(j::T) where {T <: Nemo.RingElem}
    ptr = pointer_from_objref(j)
    iptr = reinterpret(UInt, ptr) >> 24
    if !haskey(nemoNumberID, iptr)
       nemoNumberID[iptr] = live_cache(0, Array{Nemo.RingElem}(undef, 64))
    end
    val = nemoNumberID[iptr]
    val.num += 1
    push!(val.A, j)
    return reinterpret(Ptr{Cvoid}, ptr)
end

function number(j::T, j_old::T) where {T <: Nemo.RingElem}
     number_pop!(nemoNumberID, reinterpret(Ptr{Cvoid}, pointer_from_objref(j_old)))
     return number(j)
end

###############################################################################
#
#  Debug allocator
#
###############################################################################

#=
mutable struct my_counter
    added::BigInt
    deleted::BigInt
end

const line_counter_dict = Base.Dict{Ptr{Cvoid}, Int64}()
const undeletable_list = Ptr{Cvoid}[]

const pointer_compare = Any[]

const my_actual_counter = my_counter(0, 0)

function number_pop!(D::Base.Dict{UInt, live_cache}, ptr::Ptr{Cvoid})
    if haskey(line_counter_dict, ptr)
        pop!(line_counter_dict, ptr)
    else
        push!(undeletable_list, ptr)
    end
    iptr = reinterpret(UInt, ptr) >> 24
    if haskey(D, iptr)
        val = D[iptr]
        val.num -= 1
        my_actual_counter.deleted += 1
        if val.num == 0
           delete!(D, iptr)
        end
    end
end

function number(j::T, line_number) where {T <: Nemo.RingElem}
    ptr = pointer_from_objref(j)
    push!(pointer_compare, reinterpret(Ptr{Cvoid}, ptr))
    push!(pointer_compare, line_number)
    line_counter_dict[reinterpret(Ptr{Cvoid}, ptr)] = line_number
    if true
       iptr = reinterpret(UInt, ptr) >> 24
       if !haskey(nemoNumberID, iptr)
          nemoNumberID[iptr] = live_cache(0, Array{Nemo.RingElem}(undef, 64))
       end
       val = nemoNumberID[iptr]
       val.num += 1
       my_actual_counter.added += 1
       push!(val.A, j)
    end
    return reinterpret(Ptr{Cvoid},ptr)
end
=#

###############################################################################
#
#   Includes
#
###############################################################################

include("flint/fmpz.jl")
include("flint/fmpq.jl")
include("flint/fq.jl")
include("flint/fq_nmod.jl")
include("antic/nf_elem.jl")
include("nemo/Rings.jl")
include("nemo/Fields.jl")

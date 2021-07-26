# initialise a number from an mpz
function n_InitMPZ(b::BigInt, cf::coeffs_ptr)
   GC.@preserve b begin
      bb = pointer_from_objref(b)
      return n_InitMPZ_internal(bb, cf)
   end
end

# get an mpz from a number
function n_GetMPZ(s::number_ptr, r::coeffs_ptr)
   res = BigInt(1)
   GC.@preserve res begin
      resp = pointer_from_objref(res)
      n_GetMPZ_internal(resp, s, r)
   end
   return res
end

# write a number to a Singular string
function n_Write(n::number_ptr, cf::coeffs_ptr, bShortOut::Bool = false)
   d = Int(bShortOut)
   n_Write_internal(n, cf, d);
end

function n_ChineseRemainderSym(a::Vector{number_ptr}, b::Vector{number_ptr}, n::Cint, signed::Cint, cf::coeffs_ptr)
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
   omFree_internal(reinterpret(Ptr{Nothing}, m))
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
   A::Vector{Nemo.RingElem}
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

# the caller does not necessarily hold a reference to j, but julia must keep
# j alive because it is going into nemoNumberID
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

# singular requires a method to print coefficient objects as string without any
# precedence information. TODO move this to AA
function stringify_wrt_times(n)
   prec = AbstractAlgebra.prec_inf_Times
   obj = AbstractAlgebra.canonicalize(AbstractAlgebra.expressify(n))
   io = IOBuffer()
   S = AbstractAlgebra.printer(io)
   AbstractAlgebra.print_obj(S, MIME("text/plain"), obj, prec, prec)
   AbstractAlgebra.finish(S)
   return String(take!(io))
end

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


issmallinttype(::Type) = false
issmallinttype(::Type{T}) where {T<:Integer} = Base.hastypemax(T) && typemin(Int) <= typemin(T) && typemax(Int) >= typemax(T)

function number_ptr(n::Integer, c::coeffs_ptr)
   if issmallinttype(typeof(n))
      n_Init(n, c)
   else
      n_InitMPZ(BigInt(n), c)
   end
end

function number_ptr(x::Nemo.fmpz, c::coeffs_ptr)
   if Nemo._fmpz_is_small(x)
      n_Init(x.d, c)
   else
      GC.@preserve x n_InitMPZ(Nemo._as_bigint(x), c)
   end
end

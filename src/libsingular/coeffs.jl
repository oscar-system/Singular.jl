# initialise a number from an mpz
function n_InitMPZ(b::BigInt, cf::coeffs)
    bb = __mpz_struct(pointer_from_objref(b))
    return n_InitMPZ_internal(bb,cf)
end

# write a number to a Singular string
function n_Write(n::numberRef, cf::coeffs, bShortOut::Bool = false)
   d = Int(bShortOut)
   n_Write_internal(n, cf, d);
end

function n_ExtGcd(a::number, b::number, s::Ptr{number}, t::Ptr{number}, cf:: coeffs)
   sp = reinterpret(Ptr{Void},s)
   tp = reinterpret(Ptr{Void},t)
   return n_ExtGcd_internal(a, b, sp, tp, cf);
end

function n_QuotRem(a::number, b::number, p::Ptr{number}, cf::coeffs)
   pp = reinterpret(Ptr{Void}, p)
   n_QuotRem_internal(a, b, pp, cf)
end

function n_ChineseRemainderSym(a::Array{number, 1}, b::Array{number, 1}, n::Cint, signed::Cint, cf::coeffs)
   p1 = reinterpret(Ptr{Void},pointer(a))
   p2 = reinterpret(Ptr{Void},pointer(b))
   return n_ChineseRemainderSym_internal(p1, p2, n, signed, cf)
end

# create a Singular string environment
function StringSetS(m) 
   StringSetS_internal(m)
end

# omalloc free
function omFree{T}(m :: Ptr{T})
   mp = reinterpret(Ptr{Void}, m)
   omFree_internal(m)
end

###############################################################################
#
#   Setting a Ptr{number} to a number
#
###############################################################################

function setindex!(ptr::Ptr{number}, n::number)
   pp = reinterpret(Ptr{Void}, ptr)
   setindex_internal(pp, n)
end

###############################################################################
#
#   Conversion between numbers that wrap jl_value_t's and julia values
#
###############################################################################

#=

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

=#

#=
include("flint/fmpz.jl")
include("flint/fmpq.jl")
include("flint/fq.jl")
include("flint/fq_nmod.jl")
include("antic/nf_elem.jl")
include("nemo/Rings.jl")
include("nemo/Fields.jl")
=#

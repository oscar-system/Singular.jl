###############################################################################
#
#   PolyRing/spoly
#
###############################################################################

# the lookup is based on the serialized version of the ordering + nvars
# rationale: dp(0) and dp(2) look different until they meet a ring with nvars
const PolyRingID = Dict{Tuple{Union{Ring, Field}, Vector{Symbol},
                        Vector{Cint}, Int}, Nemo.MPolyRing}()

mutable struct PolyRing{T <: Nemo.RingElem} <: Nemo.MPolyRing{T}
   ptr::libSingular.ring_ptr
   refcount::Int
   base_ring::Union{Ring, Field}
   ord::sordering
   S::Vector{Symbol}

   # take ownership of the ring ptr
   function PolyRing{T}(r::libSingular.ring_ptr, R, s::Vector{Symbol}=singular_symbols(r)) where T
      @assert r.cpp_object != C_NULL
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      d = new(r, 1, R, deserialize_ordering(ord), s)
      finalizer(_PolyRing_clear_fn, d)
      return d
   end
end

function PolyRing{T}(R::Union{Ring, Field}, s::Vector{Symbol},
                     ord::sordering, cached::Bool = true,
                     degree_bound::Int = 0) where T
   bitmask = Culong(degree_bound)
   nvars = length(s)
   nvars > 0 || error("need at least one indeterminate")
   # internally in libSingular, degree_bound is set to
   deg_bound_fix = Int(libSingular.rGetExpSize(bitmask, Cint(nvars)))
   sord = serialize_ordering(nvars, ord)
   return get!(PolyRingID, (R, s, sord, deg_bound_fix)) do
      ss = rename_symbols(all_singular_symbols(R), String.(s), "x")
      v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in ss]
      r = libSingular.nCopyCoeff(R.ptr)
      ptr = libSingular.rDefault_wvhdl_helper(r, v, sord, bitmask)
      @assert deg_bound_fix == Int(libSingular.rBitmask(ptr))
      return PolyRing{T}(ptr, R, s)
   end::PolyRing{T}
end

# for various types of PolyRings, which all wrap a ring_ptr in .ptr
function _PolyRing_clear_fn(R)
   R.refcount -= 1
   if R.refcount == 0
      libSingular.rDelete(R.ptr)
   end
end

mutable struct spoly{T <: Nemo.RingElem} <: Nemo.MPolyRingElem{T}
   ptr::libSingular.poly_ptr
   parent::PolyRing{T}

   function spoly{T}(R::PolyRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

# for various types of polys, which all have a parent and a .ptr
# TODO: this is technically dirty because we can't assume that R is a valid object
function _spoly_clear_fn(p)
   R = parent(p)
   libSingular.p_Delete(p.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

function spoly{T}(R::PolyRing{T}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(0, R.ptr)
      return spoly{T}(R, p)
   end
end

function spoly{T}(R::PolyRing{T}, n::T) where T <: Nemo.RingElem
   S = parent(n)
   GC.@preserve R S n begin
      n1 = libSingular.n_Copy(n.ptr, S.ptr)
      r = libSingular.p_NSet(n1, R.ptr)
      return spoly{T}(R, r)
   end
end

# take ownership of the pointer - not for general users
function spoly{T}(R::PolyRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_NSet(n, R.ptr)
      return spoly{T}(R, p)
   end
end

function spoly{T}(R::PolyRing{T}, b::Int) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(b, R.ptr)
      return spoly{T}(R, p)
   end
end

function spoly{T}(R::PolyRing{T}, b::BigInt) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      n = libSingular.n_InitMPZ(b, S.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      return spoly{T}(R, p)
   end
end

###############################################################################
#
# Iterators
#
###############################################################################

struct SPolyCoeffs{T}
   poly::T
end

struct SPolyExponentVectors{T}
   poly::T
end

struct SPolyExponentWords{T}
   poly::T
   tmp::Vector{Int}
end

struct SPolyTerms{T}
   poly::T
end

struct SPolyMonomials{T}
   poly::T
end


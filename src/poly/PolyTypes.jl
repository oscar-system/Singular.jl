###############################################################################
#
# sordering - for fancy orderings
#
###############################################################################

# instead of many types for each block type, we have one type with possibly
# meaningless entries
struct sorder_block
   order::Singular.libSingular.rRingOrder_t
   size::Int               # per-order meaning.
   weights::Vector{Int}   # per-order meaning. empty for dp, Dp, ...
end

struct sordering
   data::Vector{sorder_block}
end

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
   base_ring::Union{Ring, Field}
   ord::sordering
   refcount::Int
   S::Vector{Symbol}

   # take ownership of a ring ptr
   function PolyRing{T}(r::libSingular.ring_ptr, R, s::Vector{Symbol}=singular_symbols(r)) where T
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      d = new(r, R, deserialize_ordering(ord), 1, s)
      finalizer(_PolyRing_clear_fn, d)
      return d
   end
end

function PolyRing{T}(R::Union{Ring, Field}, s::Vector{Symbol},
                     ord::sordering, cached::Bool = true,
                     degree_bound::Int = 0) where T
   bitmask = Culong(degree_bound)
   nvars = length(s)
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

function _PolyRing_clear_fn(R::PolyRing)
   R.refcount -= 1
   if R.refcount == 0
      libSingular.rDelete(R.ptr)
   end
end

mutable struct spoly{T <: Nemo.RingElem} <: Nemo.MPolyElem{T}
   ptr::libSingular.poly_ptr
   parent::PolyRing{T}

   function spoly{T}(R::PolyRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

function spoly{T}(R::PolyRing{T}) where T <: Nemo.RingElem
   p = libSingular.p_ISet(0, R.ptr)
   return spoly{T}(R, p)
end

function spoly{T}(R::PolyRing{T}, p::T) where T <: Nemo.RingElem
   n = libSingular.n_Copy(p.ptr, parent(p).ptr)
   r = libSingular.p_NSet(n, R.ptr)
   return spoly{T}(R, r)
end

function spoly{T}(R::PolyRing{T}, n::libSingular.number_ptr) where T <: Nemo.RingElem
   nn = libSingular.n_Copy(n, base_ring(R).ptr)
   p = libSingular.p_NSet(nn, R.ptr)
   return spoly{T}(R, p)
end

function spoly{T}(R::PolyRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   p = libSingular.p_NSet(n, R.ptr)
   return spoly{T}(R, p)
end

function spoly{T}(R::PolyRing{T}, b::Int) where T <: Nemo.RingElem
   p = libSingular.p_ISet(b, R.ptr)
   return spoly{T}(R, p)
end

function spoly{T}(R::PolyRing{T}, b::BigInt) where T <: Nemo.RingElem
   n = libSingular.n_InitMPZ(b, R.base_ring.ptr)
   p = libSingular.p_NSet(n, R.ptr)
   return spoly{T}(R, p)
end


function _spoly_clear_fn(p::spoly)
   R = parent(p)
   libSingular.p_Delete(p.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end


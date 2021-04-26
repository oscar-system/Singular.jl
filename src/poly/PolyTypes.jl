###############################################################################
#
# sordering - for fancy orderings
#
###############################################################################

struct sorder_block
   order::Singular.libSingular.rRingOrder_t
   blocksize::Int
   weights::Array{Int,1}
end

struct sordering
   data::Vector{sorder_block}
end

###############################################################################
#
#   PolyRing/spoly
#
###############################################################################

const PolyRingID = Dict{Tuple{Union{Ring, Field}, Array{Symbol, 1},
         libSingular.rRingOrder_t, libSingular.rRingOrder_t, Int}, Nemo.MPolyRing}()

mutable struct PolyRing{T <: Nemo.RingElem} <: Nemo.MPolyRing{T}
   ptr::libSingular.ring_ptr
   base_ring::Union{Ring, Field}
   ord::Union{Symbol, sordering}
   refcount::Int

   # take ownership of a ring ptr
   function PolyRing{T}(r::libSingular.ring_ptr, b, o::Symbol) where T
      d = new(r, b, o, 1)
      finalizer(_PolyRing_clear_fn, d)
      return d
   end

   # TODO these two constructors are a bit messy
   function PolyRing{T}(R::Union{Ring, Field}, s::Array{Symbol, 1},
         ord_sym::Symbol, cached::Bool = true,
         ordering::libSingular.rRingOrder_t = ringorder_dp,
         ordering2::libSingular.rRingOrder_t = ringorder_C,
         degree_bound::Int = 0) where T
      # check ordering: accept exactly one of ringorder_c, ringorder_C
      if (((ordering == ringorder_c || ordering == ringorder_C)
               && (ordering2 == ringorder_c || ordering2 == ringorder_C))
            || ((ordering != ringorder_c && ordering != ringorder_C)
               && (ordering2 != ringorder_c && ordering2 != ringorder_C)))
         error("wrong ordering")
      end
      bitmask = Culong(degree_bound)
      n_vars = Cint(length(s))
      # internally in libSingular, degree_bound is set to
      degree_bound_adjusted = Int(libSingular.rGetExpSize(bitmask, n_vars))
      if haskey(PolyRingID, (R, s, ordering, ordering2, degree_bound_adjusted))
         return PolyRingID[R, s, ordering, ordering2,
               degree_bound_adjusted]::PolyRing{T}
      else
         v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in s]
         r = libSingular.nCopyCoeff(R.ptr)

         blk0 = unsafe_wrap(Array, Ptr{Cint}(libSingular.omAlloc0(Csize_t(3*sizeof(Cint)))), 3; own=false)
         blk1 = unsafe_wrap(Array, Ptr{Cint}(libSingular.omAlloc0(Csize_t(3*sizeof(Cint)))), 3; own=false)
         if (ordering == ringorder_c || ordering == ringorder_C)
            blk0[1] = Cint(0)
            blk1[1] = Cint(0)
            blk0[2] = Cint(1)
            blk1[2] = Cint(length(v))
         else
            blk0[1] = Cint(1)
            blk1[1] = Cint(length(v))
            blk0[2] = Cint(0)
            blk1[2] = Cint(0)
         end
         ord = Array{libSingular.rRingOrder_t, 1}(undef, 3)
         ord[1] = ordering
         ord[2] = ordering2
         ord[3] = ringorder_no
         ptr = libSingular.rDefault(r, v, ord, blk0, blk1, bitmask)
         @assert degree_bound_adjusted == Int(libSingular.rBitmask(ptr))
         d = PolyRingID[R, s, ordering, ordering2, degree_bound_adjusted] =
               new(ptr, R, ord_sym, 1)
         finalizer(_PolyRing_clear_fn, d)
         return d
      end
   end

   function PolyRing{T}(R::Union{Ring, Field}, s::Array{Symbol, 1},
                        ord::sordering, cached::Bool = true,
                        degree_bound::Int = 0) where T

      bitmask = Culong(degree_bound)
      nvars = length(s)
      # internally in libSingular, degree_bound is set to
      degree_bound_adjusted = Int(libSingular.rGetExpSize(bitmask, Cint(nvars)))
         
      v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in s]
      r = libSingular.nCopyCoeff(R.ptr)
      sord = serialize_ordering(nvars, ord)
      ptr = libSingular.rDefault_wvhdl_helper(r, v, sord, bitmask)
      @assert degree_bound_adjusted == Int(libSingular.rBitmask(ptr))
      d = new(ptr, R, ord, 1)
      finalizer(_PolyRing_clear_fn, d)
      return d
   end
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
      z = new{T}(p, R)
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


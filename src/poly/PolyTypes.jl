###############################################################################
#
#   PolyRing/spoly 
#
###############################################################################

const PolyRingID = Dict{Tuple{Union{Ring, Field}, Array{Symbol, 1},
         libSingular.rRingOrder_t, libSingular.rRingOrder_t, Int}, Ring}()

mutable struct PolyRing{T <: Nemo.RingElem} <: Ring
   ptr::libSingular.ring
   base_ring::Union{Ring, Field}
   refcount::Int

   function PolyRing{T}(R::Union{Ring, Field}, s::Array{Symbol, 1},
         cached::Bool = true,
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
      n_vars = Cint(length(s));
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
               new(ptr, R, 1)
         finalizer(_PolyRing_clear_fn, d)
         return d
      end
   end
end

function _PolyRing_clear_fn(R::PolyRing)
   R.refcount -= 1
   if R.refcount == 0
      libSingular.rDelete(R.ptr)
   end
end

mutable struct spoly{T <: Nemo.RingElem} <: Nemo.RingElem
   ptr::libSingular.poly
   parent::PolyRing{T}

   function spoly{T}(R::PolyRing{T}) where T <: Nemo.RingElem
      p = libSingular.p_ISet(0, R.ptr)
	  z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
    
   function spoly{T}(R::PolyRing{T}, p::libSingular.poly) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
    
   function spoly{T}(R::PolyRing{T}, p::T) where T <: Nemo.RingElem
      n = libSingular.n_Copy(p.ptr, parent(p).ptr)
      r = libSingular.p_NSet(n, R.ptr)
      z = new(r, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
    
   function spoly{T}(R::PolyRing{T}, n::libSingular.number) where T <: Nemo.RingElem
      p = libSingular.p_NSet(n, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end

   function spoly{T}(R::PolyRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
      p = libSingular.p_NSet(n, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end

   function spoly{T}(R::PolyRing{T}, b::Int) where T <: Nemo.RingElem
      p = libSingular.p_ISet(b, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end

   function spoly{T}(R::PolyRing{T}, b::BigInt) where T <: Nemo.RingElem
      n = libSingular.n_InitMPZ(b, R.base_ring.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

function _spoly_clear_fn(p::spoly)
   R = parent(p)
   libSingular.p_Delete(p.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

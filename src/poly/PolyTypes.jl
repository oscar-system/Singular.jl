###############################################################################
#
#   PolyRing/spoly 
#
###############################################################################

const PolyRingID = Dict{Tuple{Union{Ring, Field}, Array{Symbol, 1}, libSingular.rRingOrder_t}, Ring}()

type PolyRing{T <: Nemo.RingElem} <: Ring
   ptr::libSingular.ring
   base_ring::Union{Ring, Field}
   refcount::Int

   function PolyRing{T}(R::Union{Ring, Field}, s::Array{Symbol, 1}, cached::Bool = true, 
                            ordering::libSingular.rRingOrder_t = ringorder_dp) where T
      if haskey(PolyRingID, (R, s, ordering))
         return PolyRingID[R, s, ordering]::PolyRing{T}
      else
         v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in s]
         r = libSingular.nCopyCoeff(R.ptr)
         
         ord  = unsafe_wrap(Array, 
      	   Ptr{libSingular.rRingOrder_t}(libSingular.omAlloc0(Csize_t(3*sizeof(libSingular.rRingOrder_t)))), 
      	       (Cint(3),), false)
         blk0 = unsafe_wrap(Array, Ptr{Cint}(libSingular.omAlloc0(Csize_t(3*sizeof(Cint)))), (Cint(3),), false)
         blk1 = unsafe_wrap(Array, Ptr{Cint}(libSingular.omAlloc0(Csize_t(3*sizeof(Cint)))), (Cint(3),), false)
         blk0[1] = 1
         blk1[1] = length(v)
         ord[1] = ordering
         ord[2] = ringorder_C
         ord[3] = ringorder_no
         ptr = libSingular.rDefault(r, v, ord, blk0, blk1)
         d = PolyRingID[R, s, ordering] = new(ptr, R, 1)
         finalizer(d, _PolyRing_clear_fn)
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

type spoly{T <: Nemo.RingElem} <: Nemo.RingElem
   ptr::libSingular.poly
   parent::PolyRing{T}

   function spoly{T}(R::PolyRing{T}) where T <: Nemo.RingElem
      p = libSingular.p_ISet(0, R.ptr)
	  z = new(p, R)
      R.refcount += 1
      finalizer(z, _spoly_clear_fn)
      return z
   end
    
   function spoly{T}(R::PolyRing{T}, p::libSingular.poly) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(z, _spoly_clear_fn)
      return z
   end
    
   function spoly{T}(R::PolyRing{T}, p::T) where T <: Nemo.RingElem
      n = libSingular.n_Copy(p.ptr, parent(p).ptr)
      r = libSingular.p_NSet(n, R.ptr)
      z = new(r, R)
      R.refcount += 1
      finalizer(z, _spoly_clear_fn)
      return z
   end
    
   function spoly{T}(R::PolyRing{T}, n::libSingular.number) where T <: Nemo.RingElem
      p = libSingular.p_NSet(n, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(z, _spoly_clear_fn)
      return z
   end 

   function spoly{T}(R::PolyRing{T}, b::Int) where T <: Nemo.RingElem
      p = libSingular.p_ISet(b, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(z, _spoly_clear_fn)
      return z
   end

   function spoly{T}(R::PolyRing{T}, b::BigInt) where T <: Nemo.RingElem
      n = libSingular.n_InitMPZ(b, R.base_ring.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(z, _spoly_clear_fn)
      return z
   end
end

function _spoly_clear_fn(p::spoly)
   R = parent(p)
   libSingular.p_Delete(p.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

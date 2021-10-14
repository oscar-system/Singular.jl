###############################################################################
#
#   WeylAlgebra/pweyl 
#
###############################################################################

const WeylAlgebraID = Dict{Tuple{Union{Ring, Field}, Array{Symbol, 1},
         libSingular.rRingOrder_t, libSingular.rRingOrder_t, Int},
                                                      AbstractAlgebra.NCRing}()

mutable struct WeylAlgebra{T <: Nemo.RingElem} <: AbstractAlgebra.NCRing
   ptr::libSingular.ring_ptr
   base_ring::Union{Ring, Field}
   ord::Symbol
   refcount::Int

   function WeylAlgebra{T}(R::Union{Ring, Field}, s::Array{Symbol, 1},
         ord_sym::Symbol, cached::Bool = true,
         ordering::libSingular.rRingOrder_t = ringorder_dp,
         ordering2::libSingular.rRingOrder_t = ringorder_C,
         degree_bound::Int = 0) where T

      if iszero(length(s)) || !iseven(length(s))
         error("WeylAlgebra requires an even number of symbols")
      end

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
      if haskey(WeylAlgebraID, (R, s, ordering, ordering2, degree_bound_adjusted))
         return WeylAlgebraID[R, s, ordering, ordering2,
               degree_bound_adjusted]::WeylAlgebra{T}
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
         ptr = libSingular.rWeyl(r, v, ord, blk0, blk1, bitmask)
         @assert degree_bound_adjusted == Int(libSingular.rBitmask(ptr))
         d = WeylAlgebraID[R, s, ordering, ordering2, degree_bound_adjusted] =
               new(ptr, R, ord_sym, 1)
         finalizer(_WeylAlgebra_clear_fn, d)
         return d
      end
   end
end

function (R::WeylAlgebra{T})(r::libSingular.ring_ptr) where T
    new_r = deepcopy(R)
    new_ptr = new_r.ptr
    new_r.ptr = r
    return new_r
end

function _WeylAlgebra_clear_fn(R::WeylAlgebra)
   R.refcount -= 1
   if R.refcount == 0
      libSingular.rDelete(R.ptr)
   end
end

mutable struct pweyl{T <: Nemo.RingElem} <: AbstractAlgebra.NCRingElem
   ptr::libSingular.poly_ptr
   parent::WeylAlgebra{T}

   function pweyl{T}(R::WeylAlgebra{T}) where T <: Nemo.RingElem
      p = libSingular.p_ISet(0, R.ptr)
      z = new{T}(p, R)
      R.refcount += 1
      finalizer(_pweyl_clear_fn, z)
      return z
   end
    
   function pweyl{T}(R::WeylAlgebra{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new{T}(p, R)
      R.refcount += 1
      finalizer(_pweyl_clear_fn, z)
      return z
   end
    
   function pweyl{T}(R::WeylAlgebra{T}, p::T) where T <: Nemo.RingElem
      n = libSingular.n_Copy(p.ptr, parent(p).ptr)
      r = libSingular.p_NSet(n, R.ptr)
      z = new{T}(r, R)
      R.refcount += 1
      finalizer(_pweyl_clear_fn, z)
      return z
   end
    
   function pweyl{T}(R::WeylAlgebra{T}, n::libSingular.number_ptr) where T <: Nemo.RingElem
      nn = libSingular.n_Copy(n, base_ring(R).ptr)
      p = libSingular.p_NSet(nn, R.ptr)
      z = new{T}(p, R)
      R.refcount += 1
      finalizer(_pweyl_clear_fn, z)
      return z
   end

   function pweyl{T}(R::WeylAlgebra{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
      p = libSingular.p_NSet(n, R.ptr)
      z = new(p, R)
      R.refcount += 1
      finalizer(_pweyl_clear_fn, z)
      return z
   end

   function pweyl{T}(R::WeylAlgebra{T}, b::Int) where T <: Nemo.RingElem
      p = libSingular.p_ISet(b, R.ptr)
      z = new{T}(p, R)
      R.refcount += 1
      finalizer(_pweyl_clear_fn, z)
      return z
   end

   function pweyl{T}(R::WeylAlgebra{T}, b::BigInt) where T <: Nemo.RingElem
      n = libSingular.n_InitMPZ(b, R.base_ring.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      z = new{T}(p, R)
      R.refcount += 1
      finalizer(_pweyl_clear_fn, z)
      return z
   end
end

function _pweyl_clear_fn(p::pweyl)
   R = parent(p)
   libSingular.p_Delete(p.ptr, R.ptr)
   _WeylAlgebra_clear_fn(R)
end


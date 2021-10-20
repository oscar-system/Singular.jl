###############################################################################
#
#   WeylPolyRing/sweylpoly
#
###############################################################################

const WeylPolyRingID = Dict{Tuple{Union{Ring, Field}, Vector{Symbol},
                                  Vector{Cint}, Int}, AbstractAlgebra.NCRing}()

mutable struct WeylPolyRing{T <: Nemo.RingElem} <: AbstractAlgebra.NCRing
   ptr::libSingular.ring_ptr
   refcount::Int
   base_ring::Union{Ring, Field}
   ord::sordering
   S::Vector{Symbol}

   # take ownership of a ring_ptr
   function WeylPolyRing{T}(r::libSingular.ring_ptr, R, s::Vector{Symbol}=singular_symbols(r)) where T
      @assert r.cpp_object != C_NULL
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      z = new(r, 1, R, deserialize_ordering(ord), s)
      finalizer(_PolyRing_clear_fn, z)
      return z
   end
end

function WeylPolyRing{T}(R::Union{Ring, Field}, s::Vector{Symbol},
                         ord::sordering, cached::Bool = true,
                         degree_bound::Int = 0) where T
   bitmask = Culong(degree_bound)
   nvars = length(s)
   nvars > 0 && iseven(nvars) || error("need an even number of indeterminates")
   deg_bound_fix = Int(libSingular.rGetExpSize(bitmask, Cint(nvars)))
   sord = serialize_ordering(nvars, ord)
   return get!(WeylPolyRingID, (R, s, sord, deg_bound_fix)) do
      ss = rename_symbols(all_singular_symbols(R), String.(s), "x")
      v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in ss]
      r = libSingular.nCopyCoeff(R.ptr)
      ptr = libSingular.rDefault_wvhdl_helper(r, v, sord, bitmask)
      ptr2 = libSingular.weylAlgebra(ptr)
      if libSingular.have_error()
         libSingular.rDelete(ptr2)
         error("could not construct WeylPolyRing from $ord: "*
               libSingular.get_and_clear_error())
      end
      return WeylPolyRing{T}(ptr2, R, s)
   end::WeylPolyRing{T}
end

mutable struct sweylpoly{T <: Nemo.RingElem} <: AbstractAlgebra.NCRingElem
   ptr::libSingular.poly_ptr
   parent::WeylPolyRing{T}

   function sweylpoly{T}(R::WeylPolyRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new{T}(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

function sweylpoly{T}(R::WeylPolyRing{T}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(0, R.ptr)
      return sweylpoly{T}(R, p)
   end
end

function sweylpoly{T}(R::WeylPolyRing{T}, p::T) where T <: Nemo.RingElem
   S = parent(p)
   GC.@preserve R S p begin
      n = libSingular.n_Copy(p.ptr, S.ptr)
      r = libSingular.p_NSet(n, R.ptr)
      return sweylpoly{T}(R, r)
   end
end

function sweylpoly{T}(R::WeylPolyRing{T}, n::libSingular.number_ptr) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      nn = libSingular.n_Copy(n, S.ptr)
      p = libSingular.p_NSet(nn, R.ptr)
      return sweylpoly{T}(R, p)
   end
end

function sweylpoly{T}(R::WeylPolyRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_NSet(n, R.ptr)
      return sweylpoly{T}(R, p)
   end
end

function sweylpoly{T}(R::WeylPolyRing{T}, b::Int) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(b, R.ptr)
      return sweylpoly{T}(R, p)
   end
end

function sweylpoly{T}(R::WeylPolyRing{T}, b::BigInt) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      n = libSingular.n_InitMPZ(b, S.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      return sweylpoly{T}(R, p)
   end
end


###############################################################################
#
#   ExtPolyRing/sextpoly 
#
###############################################################################

const ExtPolyRingID = Dict{Tuple{Union{Ring, Field}, Vector{Symbol},
                                 Vector{Cint}}, AbstractAlgebra.NCRing}()

mutable struct ExtPolyRing{T <: Nemo.RingElem} <: AbstractAlgebra.NCRing
   ptr::libSingular.ring_ptr
   refcount::Int
   base_ring::Union{Ring, Field}
   ord::sordering
   S::Vector{Symbol}

   # take ownership of a ring_ptr
   function ExtPolyRing{T}(r::libSingular.ring_ptr, R, s::Vector{Symbol}=singular_symbols(r)) where T <: Nemo.RingElem
      @assert r.cpp_object != C_NULL
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      z = new(r, 1, R, deserialize_ordering(ord), s)
      finalizer(_PolyRing_clear_fn, z)
      return z
   end
end

function ExtPolyRing{T}(R::Union{Ring, Field}, s::Vector{Symbol},
                        ord::sordering, cached::Bool = true) where T
   nvars = length(s)
   nvars > 1 || error("need at least two indeterminates")
   sord = serialize_ordering(nvars, ord)
   return get!(ExtPolyRingID, (R, s, sord)) do
      ss = rename_symbols(all_singular_symbols(R), String.(s), "x")
      v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in ss]
      r = libSingular.nCopyCoeff(R.ptr)
      ptr = libSingular.rDefault_wvhdl_helper(r, v, sord, Culong(0))
      ptr2 = libSingular.exteriorAlgebra(ptr)
      if libSingular.have_error()
         libSingular.rDelete(ptr2)
         error("could not construct ExtPolyRing from $ord: "*
               libSingular.get_and_clear_error())
      end
      return ExtPolyRing{T}(ptr2, R, s)
   end::ExtPolyRing{T}
end

mutable struct sextpoly{T <: Nemo.RingElem} <: AbstractAlgebra.NCRingElem
   ptr::libSingular.poly_ptr
   parent::ExtPolyRing{T}

   function sextpoly{T}(R::ExtPolyRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new{T}(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end
    
function sextpoly{T}(R::ExtPolyRing{T}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(0, R.ptr)
      return sextpoly{T}(R, p)
   end
end

function sextpoly{T}(R::ExtPolyRing{T}, p::T) where T <: Nemo.RingElem
   S = parent(p)
   GC.@preserve R S p begin
      n = libSingular.n_Copy(p.ptr, S.ptr)
      r = libSingular.p_NSet(n, R.ptr)
      return sextpoly{T}(R, r)
   end
end

function sextpoly{T}(R::ExtPolyRing{T}, n::libSingular.number_ptr) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      nn = libSingular.n_Copy(n, S.ptr)
      p = libSingular.p_NSet(nn, R.ptr)
      return sextpoly{T}(R, p)
   end
end

function sextpoly{T}(R::ExtPolyRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_NSet(n, R.ptr)
      return sextpoly{T}(R, p)
   end
end

function sextpoly{T}(R::ExtPolyRing{T}, b::Int) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(b, R.ptr)
      return sextpoly{T}(R, p)
   end
end

function sextpoly{T}(R::ExtPolyRing{T}, b::BigInt) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      n = libSingular.n_InitMPZ(b, S.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      return sextpoly{T}(R, p)
   end
end


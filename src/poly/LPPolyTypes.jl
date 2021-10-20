###############################################################################
#
#   LPPolyRing/slppoly (lp = letterplace)
#
###############################################################################

const LPPolyRingID = Dict{Tuple{Union{Ring, Field},
                                Vector{Symbol},
                                Vector{Cint},
                                Int}, AbstractAlgebra.NCRing}()

mutable struct LPPolyRing{T <: Nemo.RingElem} <: AbstractAlgebra.NCRing
   ptr::libSingular.ring_ptr
   base_ring::Union{Ring, Field}
   ord::sordering
   deg_bound::Int
   refcount::Int
   S::Vector{Symbol}

   # take ownership of a ring ptr
   function LPPolyRing{T}(r::libSingular.ring_ptr, R, deg_bound::Int,
                          s::Vector{Symbol}=singular_symbols(r)) where T
      @assert r.cpp_object != C_NULL
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      z = new(r, R, deserialize_ordering(ord), deg_bound, 1, s)
      finalizer(_PolyRing_clear_fn, z)
      return z
   end
end

function LPPolyRing{T}(R::Union{Ring, Field}, s::Vector{Symbol},
                       deg_bound::Int,
                       ord::sordering, cached::Bool = true) where T
   nvars = length(s)
   nvars > 0     || error("LPPolyRing requires at least one symbol")
   deg_bound > 0 || error("LPPolyRing requires positive degree bound")
   sord = serialize_ordering(nvars, ord)
   if cached && haskey(LPPolyRingID, (R, s, sord, deg_bound))
      return LPPolyRingID[R, s, sord, deg_bound]::LPPolyRing{T}
   else
      ss = rename_symbols(all_singular_symbols(R), String.(s), "x")
      v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in ss]
      r = libSingular.nCopyCoeff(R.ptr)
      ptr = libSingular.rDefault_wvhdl_helper(r, v, sord, Culong(1))
      ptr2 = libSingular.freeAlgebra(ptr, deg_bound, 0)
      if libSingular.have_error()
         libSingular.rDelete(ptr2)
         error("could not construct LPPolyRing from $R, $ord: "*
               libSingular.get_and_clear_error())
      end
      z = LPPolyRing{T}(ptr2, R, deg_bound, s)
      LPPolyRingID[R, s, sord, deg_bound] = z
      return z
   end
end

mutable struct slppoly{T <: Nemo.RingElem} <: Nemo.NCRingElem
   ptr::libSingular.poly_ptr
   parent::LPPolyRing{T}

   function slppoly{T}(R::LPPolyRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

function slppoly{T}(R::LPPolyRing{T}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(0, R.ptr)
      return slppoly{T}(R, p)
   end
end

function slppoly{T}(R::LPPolyRing{T}, p::T) where T <: Nemo.RingElem
   S = parent(p)
   GC.@preserve R S p begin
      n = libSingular.n_Copy(p.ptr, S.ptr)
      r = libSingular.p_NSet(n, R.ptr)
      return slppoly{T}(R, r)
   end
end

function slppoly{T}(R::LPPolyRing{T}, n::libSingular.number_ptr) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      nn = libSingular.n_Copy(n, S.ptr)
      p = libSingular.p_NSet(nn, R.ptr)
      return slppoly{T}(R, p)
   end
end

function slppoly{T}(R::LPPolyRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_NSet(n, R.ptr)
      return slppoly{T}(R, p)
   end
end

function slppoly{T}(R::LPPolyRing{T}, b::Int) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(b, R.ptr)
      return slppoly{T}(R, p)
   end
end

function slppoly{T}(R::LPPolyRing{T}, b::BigInt) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      n = libSingular.n_InitMPZ(b, S.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      return slppoly{T}(R, p)
   end
end


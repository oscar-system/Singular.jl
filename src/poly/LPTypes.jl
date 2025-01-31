###############################################################################
#
#   LPRing/slpalg (lp = letterplace)
#
###############################################################################

const LPRingID = Dict{Tuple{Union{Ring, Field},
                                Vector{Symbol},
                                Vector{Cint},
                                Int}, AbstractAlgebra.NCRing}()

mutable struct LPRing{T <: Nemo.RingElem} <: AbstractAlgebra.NCRing
   ptr::libSingular.ring_ptr
   base_ring::Union{Ring, Field}
   ord::sordering
   deg_bound::Int
   refcount::Int
   S::Vector{Symbol}

   # take ownership of the pointer - not for general users
   function LPRing{T}(r::libSingular.ring_ptr, R, deg_bound::Int,
                          s::Vector{Symbol}=singular_symbols(r)) where T
      @assert deg_bound > 0
      @assert r.cpp_object != C_NULL
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      z = new(r, R, deserialize_ordering(ord), deg_bound, 1, s)
      finalizer(_PolyRing_clear_fn, z)
      return z
   end
end

function LPRing{T}(R::Union{Ring, Field}, s::Vector{Symbol}, deg_bound::Int,
                       ord::sordering, cached::Bool = true) where T
   nvars = length(s)
   nvars > 0     || error("LPRing requires at least one symbol")
   deg_bound > 0 || error("LPRing requires positive degree bound")
   sord = serialize_ordering(nvars, ord)
   if cached && haskey(LPRingID, (R, s, sord, deg_bound))
      return LPRingID[R, s, sord, deg_bound]::LPRing{T}
   else
      ss = rename_symbols(all_singular_symbols(R), String.(s), "x")
      v = [Base.Vector{UInt8}(string(str)*"\0") for str in ss]
      r = libSingular.nCopyCoeff(R.ptr)
      ptr = GC.@preserve v libSingular.rDefault_wvhdl_helper(r, pointer.(v), sord, Culong(1))
      ptr2 = libSingular.freeAlgebra(ptr, deg_bound, 0)
      if libSingular.have_error()
         libSingular.rDelete(ptr2)
         error("could not construct LPRing from $R, $ord: "*
               libSingular.get_and_clear_error())
      end
      z = LPRing{T}(ptr2, R, deg_bound, s)
      LPRingID[R, s, sord, deg_bound] = z
      return z
   end
end

mutable struct slpalg{T <: Nemo.RingElem} <: Nemo.NCRingElem
   ptr::libSingular.poly_ptr
   parent::LPRing{T}

   # take ownership of the pointer - not for general users
   function slpalg{T}(R::LPRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

function slpalg{T}(R::LPRing{T}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(0, R.ptr)
      return slpalg{T}(R, p)
   end
end

function slpalg{T}(R::LPRing{T}, n::T) where T <: Nemo.RingElem
   S = parent(n)
   GC.@preserve R S n begin
      n1 = libSingular.n_Copy(n.ptr, S.ptr)
      r = libSingular.p_NSet(n1, R.ptr)
      return slpalg{T}(R, r)
   end
end

# take ownership of the pointer - not for general users
function slpalg{T}(R::LPRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_NSet(n, R.ptr)
      return slpalg{T}(R, p)
   end
end

function slpalg{T}(R::LPRing{T}, b::Int) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(b, R.ptr)
      return slpalg{T}(R, p)
   end
end

function slpalg{T}(R::LPRing{T}, b::BigInt) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      n = libSingular.n_InitMPZ(b, S.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      return slpalg{T}(R, p)
   end
end


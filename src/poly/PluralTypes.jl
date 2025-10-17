###############################################################################
#
#   PluralRing/spluralg
#
###############################################################################

const PluralRingID = Dict{Tuple{PolyRing, smatrix, smatrix, Vector{Symbol}},
                         AbstractAlgebra.NCRing}()

mutable struct PluralRing{T <: Nemo.RingElem} <: AbstractAlgebra.NCRing
   ptr::libSingular.ring_ptr
   refcount::Int
   base_ring::Union{Ring, Field}
   ord::sordering
   S::Vector{Symbol}

   # take ownership of the pointer - not for general users
   function PluralRing{T}(r::libSingular.ring_ptr, R, s::Vector{Symbol}=singular_symbols(r)) where T <: Nemo.RingElem
      @assert isconcretetype(T)
      @assert r.cpp_object != C_NULL
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      z = new(r, 1, R, deserialize_ordering(ord), s)
      finalizer(_PolyRing_clear_fn, z)
      return z
   end
end

function PluralRing{T}(R::PolyRing{T}, C::smatrix{spoly{T}}, D::smatrix{spoly{T}},
                      s::Vector{Symbol}, cached::Bool = true) where T <: Nemo.RingElem
   isempty(s) && error("need at least one indeterminate")
   if cached && haskey(PluralRingID, (R, C, D, s))
      return PluralRingID[R, C, D, s]::PluralRing{T}
   else
      GC.@preserve R C D begin
         ptr = libSingular.nc_CallPlural(C.ptr, D.ptr, R.ptr)
         if libSingular.have_error()
            libSingular.rDelete(ptr)
            error("could not construct G-Algebra from $R, $C, $D: "*
                   libSingular.get_and_clear_error())
         end
         @assert ptr.cpp_object != C_NULL
         return PluralRing{T}(ptr, base_ring(R), s)
      end
   end
end

mutable struct spluralg{T <: Nemo.RingElem} <: AbstractAlgebra.NCRingElem
   ptr::libSingular.poly_ptr
   parent::PluralRing{T}

   # take ownership of the pointer - not for general users
   function spluralg{T}(R::PluralRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

function spluralg{T}(R::PluralRing{T}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(0, R.ptr)
      return spluralg{T}(R, p)
   end
end

function spluralg{T}(R::PluralRing{T}, n::T) where T <: Nemo.RingElem
   S = parent(n)
   GC.@preserve R S n begin
      n1 = libSingular.n_Copy(n.ptr, S.ptr)
      r = libSingular.p_NSet(n1, R.ptr)
      return spluralg{T}(R, r)
   end
end

# take ownership of the pointer - not for general users
function spluralg{T}(R::PluralRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_NSet(n, R.ptr)
      return spluralg{T}(R, p)
   end
end

function spluralg{T}(R::PluralRing{T}, b::Int) where T <: Nemo.RingElem
   GC.@preserve R begin
      p = libSingular.p_ISet(b, R.ptr)
      return spluralg{T}(R, p)
   end
end

function spluralg{T}(R::PluralRing{T}, b::BigInt) where T <: Nemo.RingElem
   S = base_ring(R)
   GC.@preserve R S begin
      n = libSingular.n_InitMPZ(b, S.ptr)
      p = libSingular.p_NSet(n, R.ptr)
      return spluralg{T}(R, p)
   end
end


###############################################################################
#
#   GPolyRing/sgpoly
#
###############################################################################

const GPolyRingID = Dict{Tuple{PolyRing, smatrix, smatrix, Array{Symbol, 1}},
                         AbstractAlgebra.NCRing}()

mutable struct GPolyRing{T <: Nemo.RingElem} <: AbstractAlgebra.NCRing
   ptr::libSingular.ring_ptr
   base_ring::Union{Ring, Field}
   ord::sordering
   refcount::Int
   S::Vector{Symbol}

   # take ownership of a ring_ptr
   function GPolyRing{T}(r::libSingular.ring_ptr, R, s::Vector{Symbol}=singular_symbols(r)) where T <: Nemo.RingElem
      ord = Cint[]
      libSingular.rOrdering_helper(ord, r)
      d = new(r, R, deserialize_ordering(ord), 1, s)
      finalizer(_PolyRing_clear_fn, d)
      return d
   end
end

function GPolyRing{T}(R::PolyRing{T}, C::smatrix{spoly{T}}, D::smatrix{spoly{T}},
                      s::Array{Symbol, 1}, cached::Bool = true) where T <: Nemo.RingElem
   if isempty(s)
      error("GPolyRing requires at least one symbol")
   end
   if cached && haskey(GPolyRingID, (R, C, D, s))
      return GPolyRingID[R, C, D, s]::GPolyRing{T}
   else
      GC.@preserve R C D begin
         r = libSingular.nc_CallPlural(C.ptr, D.ptr, R.ptr)
         if libSingular.have_error()
            libSingular.rDelete(r)
            error("could not construct G-Algebra from $R, $C, $D: "*
                   libSingular.get_and_clear_error())
         end
         @assert r.cpp_object != C_NULL
         return GPolyRing{T}(r, base_ring(R), s)
      end
   end
end

mutable struct sgpoly{T <: Nemo.RingElem} <: AbstractAlgebra.NCRingElem
   ptr::libSingular.poly_ptr
   parent::GPolyRing{T}

   function sgpoly{T}(R::GPolyRing{T}, p::libSingular.poly_ptr) where T <: Nemo.RingElem
      z = new(p, R)
      R.refcount += 1
      finalizer(_spoly_clear_fn, z)
      return z
   end
end

function sgpoly{T}(R::GPolyRing{T}) where T <: Nemo.RingElem
   p = libSingular.p_ISet(0, R.ptr)
   return sgpoly{T}(R, p)
end

function sgpoly{T}(R::GPolyRing{T}, p::T) where T <: Nemo.RingElem
   n = libSingular.n_Copy(p.ptr, parent(p).ptr)
   r = libSingular.p_NSet(n, R.ptr)
   return sgpoly{T}(R, r)
end

function sgpoly{T}(R::GPolyRing{T}, n::libSingular.number_ptr) where T <: Nemo.RingElem
   nn = libSingular.n_Copy(n, base_ring(R).ptr)
   p = libSingular.p_NSet(nn, R.ptr)
   return sgpoly{T}(R, p)
end

function sgpoly{T}(R::GPolyRing{T}, n::Ptr{Cvoid}) where T <: Nemo.RingElem
   p = libSingular.p_NSet(n, R.ptr)
   return sgpoly{T}(R, p)
end

function sgpoly{T}(R::GPolyRing{T}, b::Int) where T <: Nemo.RingElem
   p = libSingular.p_ISet(b, R.ptr)
   return sgpoly{T}(R, p)
end

function sgpoly{T}(R::GPolyRing{T}, b::BigInt) where T <: Nemo.RingElem
   n = libSingular.n_InitMPZ(b, R.base_ring.ptr)
   p = libSingular.p_NSet(n, R.ptr)
   return sgpoly{T}(R, p)
end


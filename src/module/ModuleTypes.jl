###############################################################################
#
#   FreeMod/svector
#
###############################################################################

const FreeModID = Dict{Tuple{PolyRing, Int}, Module}()

mutable struct FreeMod{T <: Nemo.RingElem} <: Module{T}
   base_ring::PolyRing
   rank::Int

   function FreeMod{T}(R::PolyRing, r::Int) where T
      return get!(FreeModID, (R, r)) do
         new(R, r)
      end
   end
end

mutable struct svector{T <: Nemo.RingElem} <: Nemo.ModuleElem{T}
   base_ring::PolyRing
   ptr::libSingular.poly_ptr # not really a polynomial
   rank::Int

   function svector{T}(R::PolyRing, r::Int, p::libSingular.poly_ptr) where T
      z = new(R, p, r)
      R.refcount += 1
      finalizer(_svector_clear_fn, z)
      return z
   end
end

"""
    (R::PolyRing{T})(m::libSingular.poly,::Val{:vector}) where T

If R is called with a low-level poly pointer, along with
Val(:vector), it will interpret the poly pointer as a vector.
This needs to be indicated due to the fact that Singulars
vectors and polys are both stored in the poly data structure.
"""
function (R::PolyRing{T})(m::libSingular.poly_ptr, ::Val{:vector}) where T
   return svector{T}(R, 1, m)
end

function _svector_clear_fn(p::svector)
   R = p.base_ring
   libSingular.p_Delete(p.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

###############################################################################
#
#   ModuleClass/smodule
#
###############################################################################

const ModuleClassID = Dict{PolyRing, Set}()

mutable struct ModuleClass{T <: Nemo.RingElem} <: Set
   base_ring::PolyRing

   function ModuleClass{T}(R::PolyRing) where T
      return get!(ModuleClassID, R) do
         new(R)
      end
   end
end

mutable struct smodule{T <: Nemo.RingElem} <: Module{T}
   base_ring::PolyRing
   ptr::libSingular.ideal_ptr # ideal and module types are the same in Singular
   isGB::Bool

   function smodule{T}(R::PolyRing, m::libSingular.ideal_ptr) where T
      z = new(R, m, false)
      R.refcount += 1
      finalizer(_smodule_clear_fn, z)
      return z
   end
end

function smodule{T}(R::PolyRing, m::libSingular.matrix_ptr) where T
   ptr = libSingular.mp_Copy(m, R.ptr)
   ptr = libSingular.id_Matrix2Module(ptr, R.ptr)
   return smodule{T}(R, ptr)
end

function smodule{T}(R::PolyRing, vecs::svector...) where T
   n = length(vecs)
   r = vecs[1].rank;
   for i = 1:n
      @assert vecs[i].rank == r
   end
   m = libSingular.idInit(Cint(n), Cint(r))
   z = smodule{T}(R, m)
   for i = 1:n
      v = libSingular.p_Copy(vecs[i].ptr, R.ptr)
      libSingular.setindex_internal(m, v, Cint(i - 1))
   end
   return z
end

"""
    (R::PolyRing{T})(m::libSingular.ideal,::Val{:module}) where T

If R is called with a low-level ideal pointer, along with
Val(:module), it will interpret the ideal pointer as a module.
This needs to be indicated due to the fact that Singular's
modules and ideals are both stored in the ideal data structure.
"""
function (R::PolyRing{T})(m::libSingular.ideal_ptr, ::Val{:module}) where T
   return smodule{T}(R,m)
end

function _smodule_clear_fn(I::smodule)
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

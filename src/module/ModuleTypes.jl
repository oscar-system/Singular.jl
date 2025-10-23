###############################################################################
#
#   FreeMod/svector
#
###############################################################################

const FreeModID = Dict{Tuple{Nemo.NCRing, Int}, Module}()

mutable struct FreeMod{T <: Nemo.NCRingElem} <: Module{T}
   base_ring::Nemo.NCRing
   rank::Int

   function FreeMod{T}(R::Nemo.NCRing, r::Int) where T
      @assert isconcretetype(T)
      return get!(FreeModID, (R, r)) do
         new(R, r)
      end::FreeMod{T}
   end
end

mutable struct svector{T <: Nemo.NCRingElem} <: Nemo.ModuleElem{T}
   base_ring::Nemo.NCRing
   ptr::libSingular.poly_ptr # not really a polynomial
   rank::Int

   function svector{T}(R::Nemo.NCRing, r::Int, p::libSingular.poly_ptr) where T
      T === elem_type(R) || error("type mismatch")
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

function (R::PluralRing{T})(m::libSingular.poly_ptr, ::Val{:vector}) where T
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

const ModuleClassID = Dict{Nemo.NCRing, Set}()

mutable struct ModuleClass{T <: Nemo.NCRingElem} <: Set
   base_ring::Nemo.NCRing

   function ModuleClass{T}(R::Nemo.NCRing) where T
      return get!(ModuleClassID, R) do
         new(R)
      end::ModuleClass{T}
   end
end

mutable struct smodule{T <: Nemo.NCRingElem} <: Module{T}
   base_ring::Nemo.NCRing
   ptr::libSingular.ideal_ptr # ideal and module types are the same in Singular
   isGB::Bool

   # take ownership of the pointer - not for general users
   function smodule{T}(R::Nemo.NCRing, m::libSingular.ideal_ptr) where T
      T === elem_type(R) || error("type mismatch")
      z = new(R, m, false)
      R.refcount += 1
      finalizer(_smodule_clear_fn, z)
      return z
   end
end

# ownership of the pointer is NOT taken - not for general users
function smodule{T}(R::PolyRing, m::libSingular.matrix_ptr) where T
   ptr = libSingular.mp_Copy(m, R.ptr)
   ptr = libSingular.id_Matrix2Module(ptr, R.ptr)
   return smodule{T}(R, ptr)
end

function smodule{T}(R::Nemo.NCRing, vecs::svector...) where T
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
This takes ownership of the pointer and is not for general users.
"""
function (R::PolyRing)(m::libSingular.ideal_ptr, ::Val{:module})
   T = elem_type(R)
   return smodule{T}(R,m)
end

function _smodule_clear_fn(I::smodule)
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

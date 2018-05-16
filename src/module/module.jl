export smodule, ModuleClass, rank

###############################################################################
#
#   Basic manipulation 
#
###############################################################################

parent{T <: Nemo.RingElem}(a::smodule{T}) = ModuleClass{T}(a.base_ring)

base_ring(S::ModuleClass) = S.base_ring

base_ring(I::smodule) = I.base_ring

ngens(I::smodule) = I.ptr == C_NULL ? 0 : Int(libSingular.ngens(I.ptr))

#rank of the ambient space
rank(I::smodule) = Int(libSingular.rank(I.ptr))

function checkbounds(I::smodule, i::Int)
   (i > ngens(I) || i < 1) && throw(BoundsError(I, i))
end

function getindex(I::smodule{T}, i::Int) where T <: AbstractAlgebra.RingElem
   checkbounds(I, i)
   R = base_ring(I)
   p = libSingular.getindex(I.ptr, Cint(i - 1))
   return svector{T}(R, rank(I), libSingular.p_Copy(p, R.ptr))
end

iszero(p::smodule) = Bool(libSingular.idIs0(p.ptr))

function deepcopy_internal(I::smodule, dict::ObjectIdDict)
   R = base_ring(I)
   ptr = libSingular.id_Copy(I.ptr, R.ptr)
   return Module(R, ptr)
end

###############################################################################
#
#   String I/O 
#
###############################################################################

function show(io::IO, S::ModuleClass)
   print(io, "Class of Singular Modules over ")
   show(io, base_ring(S))
end

function show(io::IO, I::smodule)
   print(io, "Singular Module over ")
   show(io, base_ring(I))
   println(", with Generators:")
   n = ngens(I)
   for i = 1:n
      show(io, I[i])
      if i != n
         println(io, "")
      end
   end
end

###############################################################################
#
#   Groebner basis
#
###############################################################################

function std(I::smodule; complete_reduction::Bool=false) 
   R = base_ring(I)
   ptr = libSingular.id_Std(I.ptr, R.ptr; complete_reduction=complete_reduction)
   libSingular.idSkipZeroes(ptr)
   z = Module(R, ptr)
   z.isGB = true
   return z
end

###############################################################################
#
#   Syzygies
#
###############################################################################

function syz(M::smodule)
   R = base_ring(M)
   ptr = libSingular.id_Syzygies(M.ptr, R.ptr)
   libSingular.idSkipZeroes(ptr)
   return Module(R, ptr)
end

###############################################################################
#
#   Resolutions
#
###############################################################################

function sres{T <: Nemo.RingElem}(I::smodule{T}, max_length::Int)
   I.isGB == false && error("Not a Groebner basis ideal")
   R = base_ring(I)
   if max_length == 0
        max_length = ngens(R)
        # TODO: consider qrings
   end
   r, length, minimal = libSingular.id_sres(I.ptr, Cint(max_length + 1), R.ptr)
   for i = 1:length
      ptr = libSingular.getindex(r, Cint(i - 1))
      if ptr == C_NULL
         length = i - 1
         break
      end
     libSingular.idSkipZeroes(ptr)
   end
   return sresolution{T}(R, length, r, minimal)
end

###############################################################################
#
#   Module constructors
#
###############################################################################

function Module{T <: Nemo.RingElem}(R::PolyRing{T}, vecs::svector{T}...)
   S = elem_type(R)
   return smodule{S}(R, vecs...)
end

function Module{T <: Nemo.RingElem}(R::PolyRing{T}, id::libSingular.ideal)
   S = elem_type(R)
   return smodule{S}(R, id)
end


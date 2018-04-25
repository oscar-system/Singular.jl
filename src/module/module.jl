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

function getindex(I::smodule, i::Int)
   checkbounds(I, i)
   R = base_ring(I)
   p = libSingular.getindex(I.ptr, Cint(i - 1))
   return Singular.Vector(R, rank(I), libSingular.p_Copy(p, R.ptr))
end

iszero(p::smodule) = Bool(libSingular.idIs0(p.ptr))

function deepcopy(I::smodule)
   R = base_ring(I)
   ptr = libSingular.id_Copy(I.ptr, R.ptr)
   return SingularIdeal(R, ptr)
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
   println(io, ", with Generators:")
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
#   Resolutions
#
###############################################################################

function sres{T <: Nemo.RingElem}(I::smodule{T}, n::Int)
   I.isGB == false && error("Not a Groebner basis ideal")
   R = base_ring(I)
   len = [Cint(0)]
   r = libSingular.id_sres(I.ptr, Cint(n), pointer(len), R.ptr)
   for i = 1:n
      id = libSingular.getindex(r, Cint(i - 1))
      if id != C_NULL
         libSingular.idSkipZeroes(id)
      end
   end
   return sresolution{T}(R, n, r)
end

function lres{T <: Nemo.RingElem}(I::smodule{T}, n::Int)
   I.isGB == false && error("Not a Groebner basis ideal")
   R = base_ring(I)
   len = [Cint(n)]
   r = libSingular.id_lres(I.ptr, pointer(len), R.ptr)
   for i = 1:n
      id = libSingular.getindex(r, Cint(i - 1))
      if id != C_NULL
         libSingular.idSkipZeroes(id)
      end
   end
   return sresolution{T}(R, n, r)
end

###############################################################################
#
#   Module constructors
#
###############################################################################

function Module{T <: Nemo.RingElem}(R::PolyRing{T}, id::libSingular.ideal)
   S = elem_type(R)
   return smodule{S}(R, id)
end

# free module of rank n
function FreeModule{T <: Nemo.RingElem}(R::PolyRing{T}, n::Int)
   (n > typemax(Cint) || n < 0) && throw(DomainError())
   ptr = libSingular.id_FreeModule(Cint(n), R.ptr)
   S = elem_type(R)
   return smodule{S}(R, ptr)
end

export SingularModuleClass, SingularModule, SingularFreeModule

###############################################################################
#
#   Basic manipulation 
#
###############################################################################

parent{T <: Nemo.RingElem}(a::smodule{T}) = SingularModuleClass{T}(a.base_ring)

base_ring(S::SingularModuleClass) = S.base_ring

base_ring(I::smodule) = I.base_ring

ngens(I::smodule) = Int(libSingular.ngens(I.ptr))

function checkbounds(I::smodule, i::Int)
   (i > ngens(I) || i < 1) && throw(BoundsError(I, i))
end   

function setindex!(I::smodule, p::spoly, i::Int)
   checkbounds(I, i)
   R = base_ring(I)
   p0 = libSingular.getindex(I.ptr, Cint(i - 1))
   if p0 != C_NULL
      libSingular.p_Delete(p0, R.ptr)
   end
   p1 = libSingular.p_Copy(p.ptr, R.ptr)
   libSingular.setindex!(I.ptr, p1, Cint(i - 1))
   nothing
end

function getindex(I::smodule, i::Int)
   checkbounds(I, i)
   R = base_ring(I)
   p = libSingular.getindex(I.ptr, Cint(i - 1))
   return R(libSingular.p_Copy(p, R.ptr))
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

function show(io::IO, S::SingularModuleClass)
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
#   SingularModule constructors
#
###############################################################################

function SingularModule{T <: Nemo.RingElem}(R::SingularPolyRing{T}, id::libSingular.ideal)
   S = elem_type(R)
   return smodule{S}(R, id)
end

# free module of rank n
function SingularFreeModule{T <: Nemo.RingElem}(R::SingularPolyRing{T}, n::Int)
   (n > typemax(Cint) || n < 0) && throw(DomainError())
   ptr = libSingular.id_FreeModule(Cint(n), R.ptr)
   S = elem_type(R)
   return smodule{S}(R, ptr)
end

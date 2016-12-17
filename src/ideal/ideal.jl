export SingularIdealSet, SingularIdeal, syz

###############################################################################
#
#   Basic manipulation 
#
###############################################################################

parent{T <: Nemo.RingElem}(a::sideal{T}) = SingularIdealSet{T}(a.base_ring)

base_ring(S::SingularIdealSet) = S.base_ring

base_ring(I::sideal) = I.base_ring

ngens(I::sideal) = Int(libSingular.ngens(I.ptr))

function checkbounds(I::sideal, i::Int)
   (i > ngens(I) || i < 1) && throw(BoundsError(I, i))
end   

function setindex!(I::sideal, p::spoly, i::Int)
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

function getindex(I::sideal, i::Int)
   checkbounds(I, i)
   R = base_ring(I)
   p = libSingular.getindex(I.ptr, Cint(i - 1))
   return R(libSingular.p_Copy(p, R.ptr))
end

iszero(p::sideal) = Bool(libSingular.idIs0(p.ptr))

function deepcopy(I::sideal)
   R = base_ring(I)
   ptr = libSingular.id_Copy(I.ptr, R.ptr)
   return SingularIdeal(R, ptr)
end

###############################################################################
#
#   String I/O 
#
###############################################################################

function show(io::IO, S::SingularIdealSet)
   print(io, "Set of Singular Ideals over ")
   show(io, base_ring(S))
end

function show(io::IO, I::sideal)
   n = ngens(I)
   print(io, "Singular Ideal over ")
   show(io, base_ring(I))
   print(io, " with generators (")
   for i = 1:n
      show(io, I[i])
      if i != n
         print(io, ", ")
      end
   end
   print(io, ")")
end

###############################################################################
#
#   Groebner basis
#
###############################################################################

function std(I::sideal) 
   R = base_ring(I)
   ptr = libSingular.id_Std(I.ptr, R.ptr)
   libSingular.idSkipZeroes(ptr)
   return SingularIdeal(R, ptr)
end

###############################################################################
#
#   Groebner basis
#
###############################################################################

function syz(I::sideal) 
   R = base_ring(I)
   ptr = libSingular.id_Syzygies(I.ptr, R.ptr)
   libSingular.idSkipZeroes(ptr)
   return SingularModule(R, ptr)
end

###############################################################################
#
#   SingularIdeal constructors
#
###############################################################################

function SingularIdeal{T <: Nemo.RingElem}(R::SingularPolyRing{T}, ids::spoly{T}...)
   S = elem_type(R)
   return sideal{S}(R, ids...)
end

function SingularIdeal{T <: Nemo.RingElem}(R::SingularPolyRing{T}, id::libSingular.ideal)
   S = elem_type(R)
   return sideal{S}(R, id)
end




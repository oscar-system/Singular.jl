export SingularIdealSet, SingularIdeal, SingularMaximalIdeal, syz, lead,
       normalize!, isconstant, iszerodim

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

iszero(I::sideal) = Bool(libSingular.idIs0(I.ptr))

iszerodim(I::sideal) = Bool(libSingular.id_IsZeroDim(I.ptr, base_ring(I).ptr))

isconstant(I::sideal) = Bool(libSingular.id_IsConstant(I.ptr, base_ring(I).ptr))

function normalize!(I::sideal)
   libSingular.id_Normalize(I.ptr, base_ring(I).ptr)
   nothing
end

function deepcopy(I::sideal)
   R = base_ring(I)
   ptr = libSingular.id_Copy(I.ptr, R.ptr)
   return SingularIdeal(R, ptr)
end

function check_parent(I::sideal, J::sideal)
   base_ring(I) != base_ring(J) && error("Incompatible ideals")
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
#   Arithmetic functions
#
###############################################################################

function +(I::sideal, J::sideal)
   check_parent(I, J)
   R = base_ring(I)
   ptr = libSingular.id_Add(I.ptr, J.ptr, R.ptr)
   return SingularIdeal(R, ptr)
end

function *(I::sideal, J::sideal)
   check_parent(I, J)
   R = base_ring(I)
   ptr = libSingular.id_Mult(I.ptr, J.ptr, R.ptr)
   return SingularIdeal(R, ptr)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(I::sideal, n::Int)
   (n > typemax(Cint) || n < 0) && throw(DomainError())
   R = base_ring(I)
   ptr = libSingular.id_Power(I.ptr, Cint(n), R.ptr)
   return SingularIdeal(R, ptr)
end

###############################################################################
#
#   Leading terms
#
###############################################################################

# ideal of leading terms
function lead(I::sideal)
   R = base_ring(I)
   ptr = libSingular.id_Head(I.ptr, R.ptr)
   return SingularIdeal(R, ptr)
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

# maximal ideal in degree d
function SingularMaximalIdeal{T <: Nemo.RingElem}(R::SingularPolyRing{T}, d::Int)
   (d > typemax(Cint) || d < 0) && throw(DomainError())
   S = elem_type(R)
   ptr = libSingular.id_MaxIdeal(Cint(d), R.ptr)
   return sideal{S}(R, ptr)
end

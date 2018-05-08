export betti, minres

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(r::sresolution) = r.base_ring

base_ring(R::ResolutionSet) = R.base_ring

parent(r::sresolution) = ResolutionSet(r.base_ring)

elem_type(::Type{ResolutionSet{T}}) where T <: AbstractAlgebra.RingElem = sresolution{T}

elem_type(::ResolutionSet{T}) where T <: AbstractAlgebra.RingElem = sresolution{T}

parent_type(::Type{sresolution{T}}) where T <: AbstractAlgebra.RingElem = ResolutionSet{T}

function checkbounds(r::sresolution, i::Int)
   (i < 1 || i > r.len) && throw(BoundsError(I, i))
end

function getindex(r::sresolution, i::Int)
   checkbounds(r, i)
   R = base_ring(r)
   ptr = libSingular.getindex(r.ptr, Cint(i - 1))
   if ptr != C_NULL
      ptr = libSingular.id_Copy(ptr, R.ptr)
   end
   return Module(R, ptr)
end

length(r::sresolution) = r.len - 1

function deepcopy_internal(r::sresolution, dict::ObjectIdDict)
   R = base_ring(r)
   ptr = libSingular.res_Copy(r.ptr, r.len, R.ptr)
   return parent(r)(ptr, r.len)
end

function betti(r::sresolution)
   libSingular.syBetti(r.ptr, Cint(r.len), r.base_ring.ptr)
end

###############################################################################
#
#   Minimal resolution
#
###############################################################################

function minres(r::sresolution{T}) where T <: AbstractAlgebra.RingElem
   if r.minimal
      return r
   end
   R = base_ring(r)
   ptr = libSingular.syMinimize(r.ptr, Cint(r.len), R.ptr)
   return sresolution{T}(R, r.len, ptr, true)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::ResolutionSet)
   print(io, "Set of Singular Resolutions over ")
   show(io, R.base_ring)
end

function show(io::IO, r::sresolution)
   println(io, "Singular Resolution:")
   if r.len > 0
      ptr = libSingular.getindex(r.ptr, Cint(0))
      print(io, "R^", libSingular.rank(ptr))
   end
   for i = 1:r.len - 1
      ptr = libSingular.getindex(r.ptr, Cint(i-1))
      if ptr == C_NULL
         break
      end
      print(io, " <- R^", libSingular.ngens(ptr))
   end
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (S::ResolutionSet{T})(ptr::libSingular.resolvente, len::Int) where T <: AbstractAlgebra.RingElem
   R = base_ring(R)
   return sresolution{T}(R, len, ptr)
end


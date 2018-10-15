export ResolutionSet, sresolution, betti, minres

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(r::sresolution) = r.base_ring

base_ring(R::ResolutionSet) = R.base_ring

function parent(r::sresolution{T}) where T <: AbstractAlgebra.RingElem
   return ResolutionSet{T}(r.base_ring)
end

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
   if ptr.cpp_object != C_NULL
      ptr = libSingular.id_Copy(ptr, R.ptr)
   end
   return Module(R, ptr)
end

doc"""
    length(r::sresolution)
> Return the length of the resolution. This is what is mathematically meant by the
> length of a resolution. Over a field, this should be at most the number of variables
> in the polynomial ring.
"""
length(r::sresolution) = r.len - 1

function deepcopy_internal(r::sresolution, dict::ObjectIdDict)
   R = base_ring(r)
   ptr = libSingular.res_Copy(r.ptr, Cint(r.len), R.ptr)
   S = parent(r)
   return S(ptr, r.len)
end

###############################################################################
#
#   Betti numbers
#
###############################################################################

doc"""
    betti(r::sresolution)
> Return the Betti numbers, i.e. the ranks of the free modules in the given
> free resolution. These are returned as a Julia array of `Int`s.
"""
function betti(r::sresolution)
   libSingular.syBetti(r.ptr, Cint(r.len), r.base_ring.ptr)
end

###############################################################################
#
#   Minimal resolution
#
###############################################################################

doc"""
    minres{T <: AbstractAlgebra.RingElem}(r::sresolution{T})
> Return a minimal free resolution, given any free resolution. If the supplied
> resolution is already minimal, it may be returned without making a copy.
"""
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

function (S::ResolutionSet{T})(ptr::Ptr{Void}, len::Int) where T <: AbstractAlgebra.RingElem
   R = base_ring(S)
   return sresolution{T}(R, len, ptr)
end


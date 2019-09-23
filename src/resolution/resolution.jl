export Resolution, ResolutionSet, sresolution, betti, minres

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
   # structure has 1 more item than the mathematical length
   (i < 1 || i > length(r) + 1) && throw(BoundsError(r, i))
end

function getindex(r::sresolution, i::Int)
   checkbounds(r, i)
   R = base_ring(r)
   ptr = libSingular.getindex_internal(r.ptr, i - 1, r.minimal )
   if ptr.cpp_object != C_NULL
      ptr = libSingular.id_Copy(ptr, R.ptr)
   end
   return Module(R, ptr)
end

@doc Markdown.doc"""
    length(r::sresolution)
> Return the length of the resolution. This is what is mathematically meant by the
> length of a resolution. Over a field, this should be at most the number of variables
> in the polynomial ring.
"""
length(r::sresolution) = libSingular.get_sySize(r.ptr) - 1

function deepcopy_internal(r::sresolution, dict::IdDict)
   R = base_ring(r)
   ptr = libSingular.res_Copy(r.ptr, R.ptr)
   S = parent(r)
   return S(ptr)
end

###############################################################################
#
#   Betti numbers
#
###############################################################################

@doc Markdown.doc"""
    betti(r::sresolution)
> Return the Betti numbers, i.e. the ranks of the free modules in the given
> free resolution. These are returned as a Julia array of `Int`s. Note that the
> output of this command is useful only in the graded case.
"""
function betti(r::sresolution)
   if r.minimal
        ideal_list = libSingular.get_minimal_res(r.ptr)
   else
        ideal_list = libSingular.get_full_res(r.ptr)
   end 
   array = libSingular.syBetti_internal(ideal_list, length(r), r.base_ring.ptr)
   return unsafe_wrap(Array, array[1], array[2:3]; own=true)
end

###############################################################################
#
#   Minimal resolution
#
###############################################################################

@doc Markdown.doc"""
    minres{T <: AbstractAlgebra.RingElem}(r::sresolution{T})
> Return a minimal free resolution, given any free resolution. In the graded
> case, there exists a uniquely determined minimal resolution. If the supplied
> resolution is already minimal, it may be returned without making a copy.
"""
function minres(r::sresolution{T}) where T <: AbstractAlgebra.RingElem
   if r.minimal
      return r
   end
   R = base_ring(r)
   ptr = libSingular.syMinimize(r.ptr, R.ptr)
   return sresolution{T}(R, ptr, true)
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
   len = length(r)
   if len > 0
      ptr = libSingular.getindex_internal(r.ptr,0,r.minimal)
      if ptr.cpp_object != C_NULL
         print(io, "R^", libSingular.rank(ptr))
      end
   end
   for i = 1:len - 1
      ptr = libSingular.getindex_internal(r.ptr, i - 1, r.minimal)
      if ptr.cpp_object == C_NULL
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

function (S::ResolutionSet{T})(ptr::libSingular.syStrategy, len::Int = 0) where T <: AbstractAlgebra.RingElem
   R = base_ring(S)
   return sresolution{T}(R, ptr)
end

function (R::PolyRing{T})(ptr::libSingular.syStrategy, ::Val{:resolution}) where T <: AbstractAlgebra.RingElement
    return sresolution{T}(R, ptr, libSingular.get_minimal_res(ptr) != C_NULL )
end

###############################################################################
#
#   Resolution constructors
#
###############################################################################

@doc Markdown.doc"""
    Resolution(C::Array{smodule{T}, 1}) where T <: AbstractAlgebra.RingElem
> Create a new resolution whose maps are given by the elements of an array C of
> modules. Note that it is not checked that the maps are actually composable
> and that their pairwise composition is the zero map, that is, that the
> created resolution is a complex.
"""
function Resolution(C::Array{smodule{T}, 1}) where T <: AbstractAlgebra.RingElem
    len = size(C, 1)+1
    len > 1 || error("no module specified")
    R = base_ring(C[1])
    CC = (m -> m.ptr).(C)
    C_ptr = reinterpret(Ptr{Nothing}, pointer(CC))
    ptr = libSingular.create_SyStrategy(C_ptr, len, R.ptr)
    return sresolution{T}(R, ptr)
end


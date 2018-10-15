export smodule, ModuleClass, rank, slimgb

###############################################################################
#
#   Basic manipulation 
#
###############################################################################

parent{T <: Nemo.RingElem}(a::smodule{T}) = ModuleClass{T}(a.base_ring)

base_ring(S::ModuleClass) = S.base_ring

base_ring(I::smodule) = I.base_ring

elem_type(::ModuleClass{T}) where T <: AbstractAlgebra.RingElem = smodule{T}

elem_type(::Type{ModuleClass{T}}) where T <: AbstractAlgebra.RingElem = smodule{T}

parent_type(::Type{smodule{T}}) where T <: AbstractAlgebra.RingElem = ModuleClass{T}

doc"""
    ngens(I::smodule)
> Return the number of generators in the current representation of the module (as a list
> of vectors).
"""
ngens(I::smodule) = I.ptr == C_NULL ? 0 : Int(libSingular.ngens(I.ptr))

doc"""
    rank(I::smodule)
> Return the rank $n$ of the ambient space $R^n$ of which this module is a submodule.
"""
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

doc"""
    iszero(p::smodule)
> Return `true` if this is algebraically the zero module.
"""
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

doc"""
    std(I::smodule; complete_reduction::Bool=false)
> Compute the Groebner basis of the module $I$. If `complete_reduction` is
> set to `true`, the result is unique, up to permutation of the generators
> and multiplication by constants. If not, only the leading terms are unique
> (up to permutation of the generators and multiplication by constants, of
> course). Presently the polynomial ring used must be over a field or over
> the Singular integers.
"""
function std(I::smodule; complete_reduction::Bool=false) 
   R = base_ring(I)
   ptr = libSingular.id_Std(I.ptr, R.ptr; complete_reduction=complete_reduction)
   libSingular.idSkipZeroes(ptr)
   z = Module(R, ptr)
   z.isGB = true
   return z
end

doc"""
   slimgb(I::smodule; complete_reduction::Bool=false)
> Given a module $I$ this function computes a Groebner basis for it.
> Compared to `std`, `slimgb` uses different strategies for choosing
> a reducer.
>
> If the optional parameter `complete_reduction` is set to `true` the
> function computes a reduced Gröbner basis for $I$.
"""
function slimgb(I::smodule; complete_reduction::Bool=false)
   R = base_ring(I)
   ptr = libSingular.id_Slimgb(I.ptr, R.ptr; complete_reduction=complete_reduction)
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

doc"""
    syz(M::smodule)
> Compute the module of syzygies of the given module. This will be given as
> a set of generators in an ambient space $R^n$, where $n$ is the number of
> generators in $M$.
"""
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

doc"""
    sres{T <: Nemo.RingElem}(I::smodule{T}, max_length::Int)
> Compute a free resolution of the given module $I$ of length up to the given
> maximum length. If `max_length` is set to zero, a full length free
> resolution is computed. Each element of the resolution is itself a module.
"""
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

function Module{T <: Nemo.RingElem}(R::PolyRing{T}, vecs::svector{spoly{T}}...)
   S = elem_type(R)
   return smodule{S}(R, vecs...)
end

function Module{T <: Nemo.RingElem}(R::PolyRing{T}, id::libSingular.idealRef)
   S = elem_type(R)
   return smodule{S}(R, id)
end


export jet, minimal_generating_set, ModuleClass, rank, smodule, slimgb,
       eliminate, modulo, lift, division, divrem, prune_with_map

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(a::smodule{T}) where T <: Nemo.RingElem = ModuleClass{T}(a.base_ring)

base_ring(S::ModuleClass) = S.base_ring

base_ring(I::smodule) = I.base_ring

elem_type(::ModuleClass{T}) where T <: AbstractAlgebra.RingElem = smodule{T}

elem_type(::Type{ModuleClass{T}}) where T <: AbstractAlgebra.RingElem = smodule{T}

parent_type(::Type{smodule{T}}) where T <: AbstractAlgebra.RingElem = ModuleClass{T}


@doc raw"""
    ngens(I::smodule)

Return the number of generators in the current representation of the module (as a list
of vectors).
"""
ngens(I::smodule) = I.ptr == C_NULL ? 0 : Int(libSingular.ngens(I.ptr))

@doc raw"""
    rank(I::smodule)

Return the rank $n$ of the ambient space $R^n$ of which this module is a submodule.
"""
rank(I::smodule) = Int(GC.@preserve I libSingular.rank(I.ptr))

function checkbounds(I::smodule, i::Int)
   (i > ngens(I) || i < 1) && throw(BoundsError(I, i))
end

function getindex(I::smodule{T}, i::Int) where T <: AbstractAlgebra.RingElem
   checkbounds(I, i)
   R = base_ring(I)
   GC.@preserve I R begin
      p = libSingular.getindex(I.ptr, Cint(i - 1))
      return svector{T}(R, rank(I), libSingular.p_Copy(p, R.ptr))
   end
end

@doc raw"""
    iszero(p::smodule)

Return `true` if this is algebraically the zero module.
"""
iszero(p::smodule) = Bool(libSingular.idIs0(p.ptr))

function deepcopy_internal(I::smodule, dict::IdDict)
   R = base_ring(I)
   ptr = GC.@preserve I R libSingular.id_Copy(I.ptr, R.ptr)
   return Module(R, ptr)
end

function check_parent(I::smodule{T}, J::smodule{T}) where T <: Nemo.RingElem
   base_ring(I) != base_ring(J) && error("Incompatible modules")
end

function hash(M::smodule, h::UInt)
   v = 0x403fd5a7748e75c9%UInt
   for i in 1:ngens(M)
      v = xor(hash(M[i], h), v)
   end
   return v
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
   println(io,", with Generators:")
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
#   Leading terms
#
###############################################################################

@doc raw"""
    lead(M::smodule)

Return the module generated by the initial terms of the generators of $M$.
"""
function lead(M::smodule)
   R = base_ring(M)
   ptr = GC.@preserve M R libSingular.id_Head(M.ptr, R.ptr)
   z = Module(R, ptr)
   return z
end

###############################################################################
#
#   Groebner basis
#
###############################################################################

@doc raw"""
    std(I::smodule; complete_reduction::Bool=false)

Compute the Groebner basis of the module $I$. If `complete_reduction` is
set to `true`, the result is unique, up to permutation of the generators
and multiplication by constants. If not, only the leading terms are unique
(up to permutation of the generators and multiplication by constants, of
course). Presently the polynomial ring used must be over a field or over
the Singular integers.
"""
function std(I::smodule; complete_reduction::Bool=false)
   R = base_ring(I)
   ptr = GC.@preserve I R libSingular.id_Std(I.ptr, R.ptr, complete_reduction)
   libSingular.idSkipZeroes(ptr)
   z = Module(R, ptr)
   z.isGB = true
   return z
end

@doc raw"""
    slimgb(I::smodule; complete_reduction::Bool=false)

Given a module $I$ this function computes a Groebner basis for it.
Compared to `std`, `slimgb` uses different strategies for choosing
a reducer.

If the optional parameter `complete_reduction` is set to `true` the
function computes a reduced Gröbner basis for $I$.
"""
function slimgb(I::smodule; complete_reduction::Bool=false)
   R = base_ring(I)
   ptr = GC.@preserve I R libSingular.id_Slimgb(I.ptr, R.ptr, complete_reduction)
   libSingular.idSkipZeroes(ptr)
   z = Module(R, ptr)
   z.isGB = true
   return z
end

###############################################################################
#
#   Reduction
#
###############################################################################

@doc raw"""
    reduce(M::smodule, G::smodule; complete_reduction::Bool = true)

Return a submodule whose generators are the generators of $M$ reduced by the
submodule $G$. The submodule $G$ need not be a Groebner basis. The returned
submodule will have the same number of generators as $M$, even if they are zero.
"""
function reduce(M::smodule, G::smodule; complete_reduction::Bool = true)
   check_parent(M, G)
   R = base_ring(M)
   if complete_reduction
      ptr = GC.@preserve M G R libSingular.p_Reduce(M.ptr, G.ptr, R.ptr)
   else
      ptr = GC.@preserve M G R libSingular.p_Reduce(M.ptr, G.ptr, R.ptr,1)
   end
   return Module(R, ptr)
end

@doc raw"""
    division(I::smodule{S}, G::smodule{S}) where S

Computes a division with remainder by representing the generators of `I` in
terms of the generators of `G`. Returns a tuple (Quo, Rem, U) where
  `Matrix(I)*Matrix(U) = Matrix(G)*Matrix(Quo) + Matrix(Rem)`
where `Rem = normalform(I, std(G))`. `U` is a diagonal matrix of units differing
from the identity matrix only for local ring orderings.
"""
function division(I::smodule{S}, G::smodule{S}) where S
   check_parent(I, G)
   R = base_ring(I)
   ptr_T,ptr_Rest,ptr_U = GC.@preserve I G R libSingular.id_Lift(G.ptr, I.ptr,
                                                      true,false,true,R.ptr)
   return (smodule{S}(R,ptr_T), smodule{S}(R,ptr_Rest), smodule{S}(R,ptr_U))
end

@doc raw"""
    divrem(I::smodule{S}, G::smodule{S}; complete_reduction::Bool = false) where S <: SPolyUnion

Computes a division with remainder of the generators of `I` by
the generators of `G`. Returns a tuple (Quo, Rem, U) where
`Matrix(I)*Matrix(U) = Matrix(G)*Matrix(Quo) + Matrix(Rem)`
and `Rem = normalform(I, G)`. `U` is a diagonal matrix of units differing
from the identity matrix only for local ring orderings.
"""
function divrem(I::smodule{S}, G::smodule{S}; complete_reduction::Bool = false) where S <: SPolyUnion
   check_parent(I, G)
   R = base_ring(I)
   old_redsb=libSingular.set_option("OPT_REDSB",complete_reduction)
   old_redtail=libSingular.set_option("OPT_REDTAIL",complete_reduction)
   ptr_T,ptr_Rest,ptr_U = GC.@preserve I G R libSingular.id_Lift(G.ptr, I.ptr, true,
                                                                   false, true, R.ptr)
   libSingular.set_option("OPT_REDSB",old_redsb)
   libSingular.set_option("OPT_REDTAIL",old_redtail)
   return (smodule{S}(R,ptr_T), smodule{S}(R,ptr_Rest), smodule{S}(R,ptr_U))
end

###############################################################################
#
#   Syzygies
#
###############################################################################

@doc raw"""
    syz(M::smodule)

Compute the module of syzygies of the given module. This will be given as
a set of generators in an ambient space $R^n$, where $n$ is the number of
generators in $M$.
"""
function syz(M::smodule)
   R = base_ring(M)
   ptr = GC.@preserve M R libSingular.id_Syzygies(M.ptr, R.ptr)
   libSingular.idSkipZeroes(ptr)
   return Module(R, ptr)
end

###############################################################################
#
#   Resolutions
#
###############################################################################
@doc raw"""
    fres(id::smodule{spoly{T}}, max_length::Int, method::String="complete") where T <: Nemo.FieldElem

Compute a free resolution of the given module up to the maximum given
length. The module must be over a polynomial ring over a field, and
a Groebner basis.
The possible methods are "complete", "frame", "extended frame" and
"single module". The result is given as a resolution, whose i-th entry is
the syzygy module of the previous module, starting with the given
ideal/module.
The `max_length` can be set to $0$ if the full free resolution is required.
"""
function fres(id::smodule{spoly{T}}, max_length::Int, method::String = "complete") where T <: Nemo.FieldElem
   id.isGB || error("Not a Groebner basis")
   max_length < 0 && error("length for fres must not be negative")
   R = base_ring(id)
   if max_length == 0
        max_length = nvars(R)
        # TODO: consider qrings
   end
   if (method != "complete"
         && method != "frame"
         && method != "extended frame"
         && method != "single module")
      error("wrong optional argument for fres")
   end
   r, minimal = GC.@preserve id R libSingular.id_fres(id.ptr, Cint(max_length + 1), method, R.ptr)
   return sresolution{spoly{T}}(R, r, Bool(minimal), false)
end

@doc raw"""
    sres(I::smodule{spoly{T}}, max_length::Int) where T <: Singular.FieldElem

Compute a free resolution of the given module $I$ of length up to the given
maximum length. If `max_length` is set to zero, a full length free
resolution is computed. Each element of the resolution is itself a module.
"""
function sres(I::smodule{spoly{T}}, max_length::Int) where T <: Singular.FieldElem
   I.isGB == false && error("Not a Groebner basis ideal")
   R = base_ring(I)
   if max_length == 0
        max_length = nvars(R)
        # TODO: consider qrings
   end
   r, minimal = GC.@preserve I R libSingular.id_sres(I.ptr, Cint(max_length + 1), R.ptr)
   return sresolution{spoly{T}}(R, r, Bool(minimal))
end

@doc raw"""
    mres(id::smodule{spoly{T}}, max_length::Int) where T <: Nemo.FieldElem

Compute a minimal (free) resolution of the given module up to the maximum
given length. The module must be over a polynomial ring over a field.
The result is given as a resolution, whose i-th entry is
the syzygy module of the previous module, starting with the given module.
The `max_length` can be set to $0$ if the full free resolution is required.
"""
function mres(I::smodule{spoly{T}}, max_length::Int) where T <: Nemo.FieldElem
   R = base_ring(I)
   if max_length == 0
        max_length = nvars(R)
        # TODO: consider qrings
   end
   r, minimal = GC.@preserve I R libSingular.id_res(I.ptr, Cint(max_length + 1), 1, R.ptr)
   return sresolution{spoly{T}}(R, r, Bool(minimal), false)
end

@doc raw"""
    mres_with_map(id::smodule{spoly{T}}, max_length::Int) where T <: Nemo.FieldElem

Compute a minimal (free) resolution of the given module up to the maximum
given length. The module must be over a polynomial ring over a field.
The result is given as a resolution, whose i-th entry is
the syzygy module of the previous module, starting with the given module.
The `max_length` can be set to $0$ if the full free resolution is required.
Returns the resolution R and the transformation matrix of id to R[1].
"""
function mres_with_map(I::smodule{spoly{T}}, max_length::Int) where T <: Nemo.FieldElem
   R = base_ring(I)
   if max_length == 0
        max_length = nvars(R)
        # TODO: consider qrings
   end
   r, TT_ptr = GC.@preserve I R libSingular.id_mres_map(I.ptr, Cint(max_length + 1), R.ptr)
   return sresolution{spoly{T}}(R, r, true, false),smatrix{spoly{T}}(R,TT_ptr)
end

@doc raw"""
    prune_with_map(id::smodule{spoly{T}}) where T <: Nemo.FieldElem

Returns the module R minimally embedded in a free module such that the
corresponding factor modules are isomorphic
and the transformation matrix of id to R.
"""
function prune_with_map(I::smodule{spoly{T}}) where T <: Nemo.FieldElem
   R = base_ring(I)
   r, TT_ptr = GC.@preserve I R libSingular.id_prune_map(I.ptr, R.ptr)
   return smaodule{spoly{T}}(R, r),smatrix{spoly{T}}(R,TT_ptr)
end

@doc raw"""
    nres(id::smodule{spoly{T}}, max_length::Int) where T <: Nemo.FieldElem

Compute a minimal (free) resolution of the given module up to the maximum
given length (keeping the initial module).
The module must be over a polynomial ring over a field.
The result is given as a resolution, whose i-th entry is
the syzygy module of the previous module, starting with the given module.
The `max_length` can be set to $0$ if the full free resolution is required.
"""
function nres(I::smodule{spoly{T}}, max_length::Int) where T <: Nemo.FieldElem
   R = base_ring(I)
   if max_length == 0
        max_length = nvars(R)
        # TODO: consider qrings
   end
   r, minimal = GC.@preserve I R libSingular.id_res(I.ptr, Cint(max_length + 1), 0, R.ptr)
   return sresolution{spoly{T}}(R, r, Bool(minimal), false)
end

###############################################################################
#
#   Module constructors
#
###############################################################################

function Module(R::PolyRing{T}, vecs::svector{spoly{T}}...) where T <: Nemo.RingElem
   S = elem_type(R)
   return smodule{S}(R, vecs...)
end

function Module(R::PolyRing{T}, id::libSingular.ideal_ptr) where T <: Nemo.RingElem
   S = elem_type(R)
   return smodule{S}(R, id)
end

###############################################################################
#
#   Differential functions
#
###############################################################################

@doc raw"""
    jet(M::smodule, n::Int)

Given a module $M$ this function truncates the generators of $M$
up to degree $n$.
"""
function jet(M::smodule, n::Int)
      R = base_ring(M)
      ptr = GC.@preserve M R libSingular.id_Jet(M.ptr, Cint(n), R.ptr)
      libSingular.idSkipZeroes(ptr)
      return Module(R, ptr)
end

###############################################################################
#
#   Functions for local rings
#
###############################################################################

@doc raw"""
    minimal_generating_set(M::smodule)

Return a vector containing the minimal generators of $M$.
"""
function minimal_generating_set(M::smodule)
   R = base_ring(M)
   if has_global_ordering(R) || has_mixed_ordering(R)
      error("Ring needs local ordering.")
   end
   N = GC.@preserve M R Singular.Module(R, Singular.libSingular.idMinBase(M.ptr, R.ptr))
   return [N[i] for i in 1:ngens(N)]
end


###############################################################################
#
#   Eliminate
#
###############################################################################

@doc raw"""
    eliminate(M::smodule, polys::spoly...)

Given a module `M` and a list of polynomials which are variables, construct the
intersection of `M` with the free module where those variables have been eliminated.
"""
function eliminate(M::smodule, polys::spoly...)
   R = base_ring(M)
   p = one(R)
   for i = 1:length(polys)
      !is_gen(polys[i]) && error("Not a variable")
      parent(polys[i]) != R && error("Incompatible base rings")
      p *= polys[i]
   end
   ptr = GC.@preserve M p R libSingular.id_Eliminate(M.ptr, p.ptr, R.ptr)
   return Module(R, ptr)
end

###############################################################################
#
#   Lift
#
###############################################################################

@doc raw"""
    lift(M::smodule{T}, SM::smodule{T}) where T

Represents the generators of `SM` in terms of the generators of `M`.
If `SM` is in `M`, `rest` is the null module, otherwise `rest = reduce(SM, std(M))`.
Returns `(result, rest)` with
for global orderings:
    `Matrix(SM) - Matrix(rest) = Matrix(M)*Matrix(result)`
for non-global orderings:
    `Matrix(SM)*U - Matrix(rest) = Matrix(M)*Matrix(result)`
where `U` is some diagonal matrix of units. To compute this `U`,
see `lift(M::smodule, SM::smodule, goodShape::Bool, isSB::Bool, divide::Bool)`.
"""
function lift(M::smodule{T}, SM::smodule{T}) where T
   R = base_ring(M)
   R == base_ring(SM) || error("base rings must match")
   ptr, rest_ptr = GC.@preserve M SM R libSingular.id_Lift(M.ptr, SM.ptr, R.ptr)
   return Module(R, ptr), Module(R,rest_ptr)
end

@doc raw"""
    lift(M::smodule{T}, SM::smodule{T}, goodShape::Bool, isSB::Bool, divide::Bool) where T

Represents the generators of `SM` in terms of the generators of `M`.
Returns `(result, rest, U)` with
    `Matrix(SM)*U - Matrix(rest) = Matrix(M)*Matrix(result)`
If `SM` is in `M`, then `rest` is the null module. Otherwise, `rest = SM` if
`!divide`, and `rest = normalform(SM, std(M))` if `divide`.
`U` is a diagonal matrix of units, differing from the identity matrix only for
local ring orderings.

There are three boolean options.
`goodShape`: maximal non-zero index in generators of `SM` <= that of `M`, which
should be come from a rank check `rank(SM)==rank(M)`.
`isSB`: generators of `M` form a Groebner basis.
`divide`: allow `SM` not to be a submodule of `M`.
"""
function lift(M::smodule{T}, SM::smodule{T},
                            goodShape::Bool, isSB::Bool, divide::Bool) where T
   R = base_ring(M)
   R == base_ring(SM) || error("base rings must match")
   res, rest, U = GC.@preserve M SM R libSingular.id_Lift(M.ptr, SM.ptr,
                                                goodShape, isSB, divide, R.ptr)
   return (smodule{T}(R, res), smodule{T}(R, rest), smatrix{T}(R, U))
end

###############################################################################
#
#   LiftStd
#
###############################################################################

@doc raw"""
    lift_std_syz(M::smodule)

Computes the Groebner base `G` of `M`, the transformation matrix `T` and the syzygies of M.
Returns a tuple `(G,T,S)` satisfying `(Matrix(G) = Matrix(M) * T, 0=Matrix(M)*Matrix(S))`.
"""
function lift_std_syz(M::smodule; complete_reduction::Bool = false)
   R = base_ring(M)
   ptr,T_ptr,S_ptr = GC.@preserve M R libSingular.id_LiftStdSyz(M.ptr, R.ptr, complete_reduction)
   return Module(R, ptr), smatrix{elem_type(R)}(R, T_ptr), Module(R,S_ptr)
end

@doc raw"""
    lift_std(M::smodule)

Computes the Groebner base `G` of `M` and the transformation matrix `T` such that
`(Matrix(G) = Matrix(M) * T)`.
"""
function lift_std(M::smodule; complete_reduction::Bool = false)
   R = base_ring(M)
   ptr,T_ptr = GC.@preserve M R libSingular.id_LiftStd(M.ptr, R.ptr, complete_reduction)
   return Module(R, ptr), smatrix{elem_type(R)}(R, T_ptr)
end

###############################################################################
#
#   Modulo
#
###############################################################################

@doc raw"""
    modulo(A::smodule, B:smodule)

Represents  A/(A intersect B) (isomorphic to (A+B)/B)
"""
function modulo(A::smodule, B::smodule)
   R = base_ring(A)
   ptr = GC.@preserve A B R libSingular.id_Modulo(A.ptr, B.ptr, R.ptr)
   return Module(R, ptr)
end


function vdim(I::smodule)
   I.isGB || error("Not a Groebner basis")
   R = base_ring(I)
   GC.@preserve I R return Int(libSingular.id_vdim(I.ptr, R.ptr))
end

###############################################################################
#
#   Hilbert series
#
###############################################################################

function hilbert_series(M::smodule{spoly{T}}, w::Vector{<:Integer}, shifts::Vector{<:Integer}) where T <: Nemo.FieldElem
  M.isGB || error("Not a Groebner basis")
  R = base_ring(M)
  length(w) == nvars(R) || error("wrong number of weights")
  length(shifts) == rank(M) || error("wrong number of weights")
  all(x -> x>0, w) || error("weights must be positive")
  w = convert(Vector{Int32}, w)
  shifts = convert(Vector{Int32}, shifts)
  z = Vector{Int32}()
  GC.@preserve M R libSingular.scHilbWeighted(M.ptr, R.ptr, w, shifts, z)
  return z
end

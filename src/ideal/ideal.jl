export sideal, IdealSet, syz, lead, normalize!, isconstant, iszerodim, fres,
       sres, lres, intersection, quotient, reduce, eliminate, kernel, equal,
       contains

###############################################################################
#
#   Basic manipulation 
#
###############################################################################

parent{T <: Nemo.RingElem}(a::sideal{T}) = IdealSet{T}(a.base_ring)

base_ring(S::IdealSet) = S.base_ring

base_ring(I::sideal) = I.base_ring

elem_type(::Type{IdealSet{spoly{T}}}) where T <: Nemo.RingElem = sideal{spoly{T}}

elem_type(::IdealSet{spoly{T}}) where T <: Nemo.RingElem = sideal{spoly{T}}

parent_type(::Type{sideal{spoly{T}}}) where T <: Nemo.RingElem = IdealSet{spoly{T}}

ngens(I::sideal) = Int(libSingular.ngens(I.ptr))

function checkbounds(I::sideal, i::Int)
   (i > ngens(I) || i < 1) && throw(BoundsError(I, i))
end   

function setindex!{T <: Nemo.RingElem}(I::sideal{spoly{T}}, p::spoly{T}, i::Int)
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
   return Ideal(R, ptr)
end

function check_parent{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})
   base_ring(I) != base_ring(J) && error("Incompatible ideals")
end

###############################################################################
#
#   String I/O 
#
###############################################################################

function show(io::IO, S::IdealSet)
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

function +{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})
   check_parent(I, J)
   R = base_ring(I)
   ptr = libSingular.id_Add(I.ptr, J.ptr, R.ptr)
   return Ideal(R, ptr)
end

function *{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})
   check_parent(I, J)
   R = base_ring(I)
   ptr = libSingular.id_Mult(I.ptr, J.ptr, R.ptr)
   return Ideal(R, ptr)
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
   return Ideal(R, ptr)
end

###############################################################################
#
#   Containment
#
###############################################################################

doc"""
    contains{T <: AbstractAlgebra.RingElem}(I::sideal, J::sideal)
> Returns `true` if the ideal $I$ contains the ideal $J$. This will be
> expensive if $I$ is not a Groebner ideal, since its standard form must be
> computed.
"""
function contains(I::sideal{T}, J::sideal{T}) where T <: AbstractAlgebra.RingElem
   if !I.isGB
      I = std(I)
   end
   return iszero(reduce(J, I))
end

###############################################################################
#
#   Comparison
#
###############################################################################

doc"""
    isequal{T <: AbstractAlgebra.RingElem}(I1::sideal{T}, I2::sideal{T})
> Return `true` if the given ideals have the same generators in the same order. Note
> that two arithmetically equal ideals with different generators will return `false`.
"""
function isequal(I1::sideal{T}, I2::sideal{T}) where T <: AbstractAlgebra.RingElem
   check_parent(I1, I2)
   if ngens(I1) != ngens(I2)
      return false
   end
   R = base_ring(I1)
   return Bool(libSingular.id_IsEqual(I1.ptr, I2.ptr, R.ptr))
end

doc"""
    equal(I1::sideal{T}, I2::sideal{T}) where T <: AbstractAlgebra.RingElem
> Return `true` if the two ideals are contained in each other, i.e. are the same
> ideal mathematically. This function should be called only as a last
> resort; it is exceptionally expensive to test equality of ideals! Do not
> define `==` as an alias for this function!
"""
function equal(I1::sideal{T}, I2::sideal{T}) where T <: AbstractAlgebra.RingElem
   check_parent(I1, I2)
   return contains(I1, I2) && contains(I2, I1)
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
   return Ideal(R, ptr)
end

###############################################################################
#
#   Intersection
#
###############################################################################

function intersection{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})
   check_parent(I, J)
   R = base_ring(I)
   ptr = libSingular.id_Intersection(I.ptr, J.ptr, R.ptr)
   return Ideal(R, ptr)
end

###############################################################################
#
#   Quotient
#
###############################################################################

function quotient{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})
   check_parent(I, J)
   R = base_ring(I)
   ptr = libSingular.id_Quotient(I.ptr, J.ptr, I.isGB, R.ptr)
   return Ideal(R, ptr)
end

###############################################################################
#
#   Groebner basis
#
###############################################################################

function std(I::sideal; complete_reduction::Bool=false) 
   R = base_ring(I)
   ptr = libSingular.id_Std(I.ptr, R.ptr; complete_reduction=complete_reduction)
   libSingular.idSkipZeroes(ptr)
   z = Ideal(R, ptr)
   z.isGB = true
   return z
end

function satstd(I::sideal, J::sideal)
   check_parent(I, J)
   R = base_ring(I)
   ptr = libSingular.id_Satstd(I.ptr, J.ptr, R.ptr)
   libSingular.idSkipZeroes(ptr)
   z = Ideal(R, ptr)
   z.isGB = true
   return z
end

###############################################################################
#
#   Reduction
#
###############################################################################

function reduce(I::sideal, G::sideal)
   check_parent(I, G)
   R = base_ring(I)
   !G.isGB && error("Not a Groebner basis")
   ptr = libSingular.p_Reduce(I.ptr, G.ptr, R.ptr)
   return Ideal(R, ptr)
end

function reduce(p::spoly, G::sideal)
   R = base_ring(G)
   par = parent(p)
   R != par && error("Incompatible base rings")
   !G.isGB && error("Not a Groebner basis")
   ptr = libSingular.p_Reduce(p.ptr, G.ptr, R.ptr)
   return par(ptr)
end

###############################################################################
#
#   Eliminate
#
###############################################################################

function eliminate(I::sideal, p::spoly)
   R = base_ring(I)
   parent(p) != R && error("Incompatible base rings")
   ptr = libSingular.id_Eliminate(I.ptr, p.ptr, R.ptr)
   return Ideal(R, ptr)
end

#=
The kernel of the map \phi defined as follows:
Let v_1, ..., v_s be the variables in the polynomial ring 'source'. Then
\phi(v_i) := map[i].
This is internally computed via elimination.
=#
function kernel(source::PolyRing, map::sideal)
   # TODO: check for quotient rings and/or local (or mixed) orderings, see
   #       jjPREIMAGE() in the Singular interpreter
   target = base_ring(map)
   zero_ideal = Ideal(target, )
   ptr = libSingular.maGetPreimage(target.ptr, map.ptr, zero_ideal.ptr,
         source.ptr)
   return Ideal(source, ptr)
end

###############################################################################
#
#   Syzygies
#
###############################################################################

function syz(I::sideal) 
   R = base_ring(I)
   ptr = libSingular.id_Syzygies(I.ptr, R.ptr)
   libSingular.idSkipZeroes(ptr)
   return Module(R, ptr)
end

###############################################################################
#
#   Resolutions
#
###############################################################################

function fres{T <: Nemo.RingElem}(id::sideal{T}, max_length::Int,
      method::String="complete")
   id.isGB == false && error("ideal is not a standard basis")
   max_length < 0 && error("length for fres must not be negative")
   R = base_ring(id)
   if max_length == 0
        max_length = ngens(R)+1
        # TODO: consider qrings
   end
   if (method != "complete"
         && method != "frame"
         && method != "extended frame"
         && method != "single module")
      error("wrong optional argument for fres")
   end
   r, length = libSingular.id_fres(id.ptr, Cint(max_length), method, R.ptr)
   return sresolution{T}(R, length, r)
end

function sres{T <: Nemo.RingElem}(I::sideal{T}, n::Int)
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

function lres{T <: Nemo.RingElem}(I::sideal{T}, n::Int)
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
#   Ideal constructors
#
###############################################################################

function Ideal{T <: Nemo.RingElem}(R::PolyRing{T}, ids::spoly{T}...)
   S = elem_type(R)
   return sideal{S}(R, ids...)
end

function Ideal{T <: Nemo.RingElem}(R::PolyRing{T}, ids::Array{spoly{T}, 1})
   S = elem_type(R)
   return sideal{S}(R, ids...)
end


function Ideal{T <: Nemo.RingElem}(R::PolyRing{T}, id::libSingular.ideal)
   S = elem_type(R)
   return sideal{S}(R, id)
end

# maximal ideal in degree d
function MaximalIdeal{T <: Nemo.RingElem}(R::PolyRing{T}, d::Int)
   (d > typemax(Cint) || d < 0) && throw(DomainError())
   S = elem_type(R)
   ptr = libSingular.id_MaxIdeal(Cint(d), R.ptr)
   return sideal{S}(R, ptr)
end

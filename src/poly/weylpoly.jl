export sweylpoly, WeylPolyRing, WeylAlgebra

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::WeylPolyRing{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::sweylpoly) = base_ring(parent(p))

elem_type(::Type{WeylPolyRing{T}}) where T <: Nemo.RingElem = sweylpoly{T}

parent_type(::Type{sweylpoly{T}}) where T <: Nemo.RingElem = WeylPolyRing{T}

@doc Markdown.doc"""
    degree_bound(R::WeylPolyRing)

Return the internal degree bound in each variable, enforced by Singular. This is the
largest positive value any degree can have before an overflow will occur. This
internal bound may be higher than the bound requested by the user via the
`degree_bound` parameter of the `WeylAlgebra` constructor.
"""
function degree_bound(R::WeylPolyRing)
   return Int(libSingular.rBitmask(R.ptr))
end


###############################################################################
#
#   String I/O
#
###############################################################################          

function show(io::IO, R::WeylPolyRing)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
     print(io, "Singular Weyl Algebra Quotient Ring ", s)
   else
     print(io, "Singular Weyl Algebra ", s)
   end
end

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{sweylpoly{T}}, ::Type{sweylpoly{T}}) where T <: Nemo.RingElem = sweylpoly{T}

function promote_rule(::Type{sweylpoly{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? sweylpoly{T} : Union{}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::WeylPolyRing)()
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R)
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::Int)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, ptr)
   z.parent = R
   return z
end

function (R::WeylPolyRing)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   z = sweylpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::WeylPolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Weyl algebra")
   z = sweylpoly{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::WeylPolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::WeylPolyRing)(p::sweylpoly)
   parent(p) != R && error("Unable to coerce")
   return p
end

function (R::WeylPolyRing)(n::libSingular.number_ptr)
    return R.base_ring(n)
end

###############################################################################
#
#   WeylAlgebra constructor
#
###############################################################################

function _WeylAlgebra(R, s::Vector{String}, ordering, ordering2, cached, degree_bound)
   sord = get_fancy_ordering(ordering, ordering2)
   z = WeylPolyRing{elem_type(R)}(R, Symbol.(s), sord, cached, degree_bound)
   return (z, gens(z))
end

function WeylAlgebra(R::Union{Ring, Field}, s::Vector{String};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   s = vcat(s, ["d"*sym for sym in s])
   return _WeylAlgebra(R, s, ordering, ordering2, cached, degree_bound)
end

function WeylAlgebra(R::Union{Ring, Field}, s::Matrix{String};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   s = vcat(view(s, 1, :), view(s, 2, :))
   return _WeylAlgebra(R, s, ordering, ordering2, cached, degree_bound)
end

function WeylAlgebra(R::Nemo.Ring, s::Vector{String};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   R = CoefficientRing(R)
   s = vcat(s, ["d"*sym for sym in s])
   return _WeylAlgebra(R, s, ordering, ordering2, cached, degree_bound)
end

function WeylAlgebra(R::Nemo.Ring, s2::Matrix{String};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   R = CoefficientRing(R)
   s = vcat(view(s, 1, :), view(s, 2, :))
   return _WeylAlgebra(R, s, ordering, ordering2, cached, degree_bound)
end

macro WeylAlgebra(R, s, n, o)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = WeylAlgebra($R, $v0; ordering=$o))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

macro WeylAlgebra(R, s, n)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = WeylAlgebra($R, $v0))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

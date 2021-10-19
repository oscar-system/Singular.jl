export sextpoly, ExtPolyRing, ExteriorAlgebra

###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(R::ExtPolyRing{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::sextpoly) = base_ring(parent(p))

elem_type(::Type{ExtPolyRing{T}}) where T <: Nemo.RingElem = sextpoly{T}

parent_type(::Type{sextpoly{T}}) where T <: Nemo.RingElem = ExtPolyRing{T}

nvars(R::ExtPolyRing) = Int(libSingular.rVar(R.ptr))

has_global_ordering(R::ExtPolyRing) = Bool(libSingular.rHasGlobalOrdering(R.ptr))

has_mixed_ordering(R::ExtPolyRing) = Bool(libSingular.rHasMixedOrdering(R.ptr))

function has_local_ordering(R::ExtPolyRing)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

isquotient_ring(R::ExtPolyRing) = Bool(Singular.libSingular.rIsQuotientRing(R.ptr))

characteristic(R::ExtPolyRing) = Int(libSingular.rChar(R.ptr))

function gens(R::ExtPolyRing)
   n = nvars(R)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::ExtPolyRing, i::Int)
   return R(libSingular.rGetVar(Cint(i), R.ptr))
end

function symbols(R::ExtPolyRing)
   io = IOBuffer();
   symbols = Array{Symbol,1}(undef, 0)
   for g in gens(R)
      show(io,g)
      push!(symbols, Symbol(String(take!(io))))
   end
   return symbols
end

ordering(R::ExtPolyRing) = R.ord

function degree_bound(R::ExtPolyRing)
   return Int(libSingular.rBitmask(R.ptr))
end

zero(R::ExtPolyRing) = R()

one(R::ExtPolyRing) = R(1)


function Base.hash(p::sextpoly{T}, h::UInt) where T <: Nemo.RingElem
   v = 0xaf708b07f940b4d2%UInt
   v = xor(hash(collect(exponent_vectors(p)), h), v)
   for c in coefficients(p)
      v = xor(hash(c, h), v)
      v = (v << 1) | (v >> (sizeof(Int)*8 - 1))
   end
   return v
end

###############################################################################
#
#   String I/O
#
###############################################################################          
function show(io::IO, R::ExtPolyRing)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
     print(io, "Singular Exterior Algebra Quotient Ring ", s)
   else
     print(io, "Singular Exterior Algebra ", s)
   end
end

show_minus_one(::Type{sextpoly{T}}) where T <: Nemo.RingElem = show_minus_one(T)

needs_parentheses(x::sextpoly) = length(x) > 1

isnegative(x::sextpoly) = isconstant(x) && !iszero(x) && isnegative(coeff(x, 0))

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{sextpoly{T}}, ::Type{sextpoly{T}}) where T <: Nemo.RingElem = sextpoly{T}

function promote_rule(::Type{sextpoly{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? sextpoly{T} : Union{}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::ExtPolyRing)()
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R)
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::Int)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, ptr)
   z.parent = R
   return z
end

function (R::ExtPolyRing)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   z = sextpoly{T}(R, n)
   z.parent = R
   return z
end

function (R::ExtPolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Exterior algebra")
   z = sextpoly{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::ExtPolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::ExtPolyRing)(p::sextpoly)
   parent(p) != R && error("Unable to coerce")
   return p
end

function(R::ExtPolyRing)(n::libSingular.number_ptr)
    return R.base_ring(n)
end

###############################################################################
#
#   ExtPolyRing constructor
#
###############################################################################

function ExteriorAlgebra(R::Union{Ring, Field}, s::Array{String, 1};
                         cached::Bool = true,
                         ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min)
   U = [Symbol(v) for v in s]
   T = elem_type(R)
   z = ExtPolyRing{T}(R, U, ordering, cached, sym2ringorder[ordering],
                                              sym2ringorder[ordering2])
   return z, gens(z)
end

function ExteriorAlgebra(R::Nemo.Ring, s::Array{String, 1};
                         cached::Bool = true,
                         ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min)
   S = CoefficientRing(R)
   U = [Symbol(v) for v in s]
   T = elem_type(S)
   z = ExtPolyRing{T}(S, U, ordering, cached, sym2ringorder[ordering],
                                              sym2ringorder[ordering2])
   return z, gens(z)
end

macro ExteriorAlgebra(R, s, n, o)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = ExtPolyRing($R, $v0; ordering=$o))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

macro ExteriorAlgebra(R, s, n)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = ExtPolyRing($R, $v0))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

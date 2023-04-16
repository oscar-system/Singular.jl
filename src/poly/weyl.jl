export WeylAlgebra

const WeylAlgebraID = Dict{Tuple{Union{Ring, Field}, Vector{Symbol},
                                  Vector{Cint}, Int}, AbstractAlgebra.NCRing}()

function _WeylAlgebra(R, s::AbstractVector{<:VarName}, ordering, ordering2, cached, degree_bound)
   return _WeylAlgebra(R, map(Symbol, s), ordering, ordering2, cached, degree_bound)
end
function _WeylAlgebra(R, s::Vector{Symbol}, ordering, ordering2, cached, degree_bound)
   bitmask = Culong(degree_bound)
   nvars = length(s)
   nvars > 0 && iseven(nvars) || error("need an even number of indeterminates")
   ord = get_fancy_ordering(ordering, ordering2)
   deg_bound_fix = Int(libSingular.rGetExpSize(bitmask, Cint(nvars)))
   sord = serialize_ordering(nvars, ord)
   T = elem_type(R)
   z = get!(WeylAlgebraID, (R, s, sord, deg_bound_fix)) do
            ss = rename_symbols(all_singular_symbols(R), String.(s), "x")
            v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in ss]
            r = libSingular.nCopyCoeff(R.ptr)
            ptr = libSingular.rDefault_wvhdl_helper(r, v, sord, bitmask)
            ptr2 = libSingular.weylAlgebra(ptr)
            if libSingular.have_error()
               libSingular.rDelete(ptr2)
               error("could not construct WeylAlgebra from $ord: "*
                     libSingular.get_and_clear_error())
            end
            return PluralRing{T}(ptr2, R, s)
         end::PluralRing{T}
   return (z, gens(z))
end

function WeylAlgebra(R::Union{Ring, Field}, s::AbstractVector{<:VarName};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   s = vcat(map(Symbol, s), [Symbol('d', sym) for sym in s])
   return _WeylAlgebra(R, s, ordering, ordering2, cached, degree_bound)
end

function WeylAlgebra(R::Union{Ring, Field}, s::AbstractMatrix{<:VarName};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   size(s)[1] == 2 || throw(ArgumentError("s must be either a vector or a matrix with two rows"))
   s = vcat(view(s, 1, :), view(s, 2, :))
   return _WeylAlgebra(R, map(Symbol, s), ordering, ordering2, cached, degree_bound)
end

function WeylAlgebra(R::Nemo.Ring, s::AbstractVector{<:VarName};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   s = vcat(map(Symbol, s), [Symbol('d', sym) for sym in s])
   RR = CoefficientRing(R)
   return _WeylAlgebra(RR, map(Symbol, s), ordering, ordering2, cached, degree_bound)
end

function WeylAlgebra(R::Nemo.Ring, s::AbstractMatrix{<:VarName};
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true, degree_bound::Int = 0)
   size(s)[1] == 2 || throw(ArgumentError("s must be either a vector or a matrix with two rows"))
   s = vcat(view(s, 1, :), view(s, 2, :))
   RR = CoefficientRing(R)
   return _WeylAlgebra(RR, map(Symbol, s), ordering, ordering2, cached, degree_bound)
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


###############################################################################
#
#   Fancy orderings
#
###############################################################################

function AbstractAlgebra.expressify(a::sordering; context = nothing)
   prod = Expr(:call, :cdot)
   for i in a.data
      if i.order == ringorder_lp
         this = Expr(:call, :ordering_lp, i.size)
      elseif i.order == ringorder_rp
         this = Expr(:call, :ordering_rp, i.size)
      elseif i.order == ringorder_dp
         this = Expr(:call, :ordering_dp, i.size)
      elseif i.order == ringorder_Dp
         this = Expr(:call, :ordering_Dp, i.size)
      elseif i.order == ringorder_wp
         this = Expr(:call, :ordering_wp, string(i.weights))
      elseif i.order == ringorder_Wp
         this = Expr(:call, :ordering_Wp, string(i.weights))
      elseif i.order == ringorder_ls
         this = Expr(:call, :ordering_ls, i.size)
      elseif i.order == ringorder_rs
         this = Expr(:call, :ordering_rs, i.size)
      elseif i.order == ringorder_ds
         this = Expr(:call, :ordering_ds, i.size)
      elseif i.order == ringorder_Ds
         this = Expr(:call, :ordering_Ds, i.size)
      elseif i.order == ringorder_ws
         this = Expr(:call, :ordering_ws, string(i.weights))
      elseif i.order == ringorder_Ws
         this = Expr(:call, :ordering_Ws, string(i.weights))
      elseif i.order == ringorder_a
         this = Expr(:call, :ordering_a, string(i.weights))
      elseif i.order == ringorder_M
         this = Expr(:call, :ordering_M, string(transpose(reshape(i.weights, (i.size, i.size)))))
      elseif i.order == ringorder_C
         this = Expr(:call, :ordering_C)
      elseif i.order == ringorder_c
         this = Expr(:call, :ordering_c)
      elseif i.order == ringorder_S
         this = Expr(:call, :ordering_S)
      elseif i.order == ringorder_s
         this = Expr(:call, :ordering_s, i.size)
      else
         this = Expr(:call, :ordering_unknown)
      end
      push!(prod.args, this)
   end
   return prod
end

function Base.show(io::IO, mi::MIME"text/plain", a::sordering)
   Singular.AbstractAlgebra.show_via_expressify(io, mi, a)
end

function Base.show(io::IO, a::sordering)
   Singular.AbstractAlgebra.show_via_expressify(io, a)
end

function _is_basic_ordering(t::libSingular.rRingOrder_t)
    return t == ringorder_lp || t == ringorder_ls ||
           t == ringorder_rp || t == ringorder_rs ||
           t == ringorder_dp || t == ringorder_ds ||
           t == ringorder_Dp || t == ringorder_Ds
end

function _is_weighted_ordering(t::libSingular.rRingOrder_t)
    return t == ringorder_wp || t == ringorder_ws ||
           t == ringorder_Wp || t == ringorder_Ws
end

function _basic_ordering(t::libSingular.rRingOrder_t, size::Int)
   size >= 0 || throw(ArgumentError("block size must be nonnegative"))
   return sordering([sorder_block(t, size, Int[])])
end

function _global_weighted_ordering(t::libSingular.rRingOrder_t, v::Vector{Int})
   len = length(v)
   len > 0 || throw(ArgumentError("weight vector must be non-empty"))
   all(x->x>0, v) || throw(ArgumentError("all weights must be positive"))
   return sordering([sorder_block(t, len, v)])
end

function _local_weighted_ordering(t::libSingular.rRingOrder_t, v::Vector{Int})
   len = length(v)
   len > 0 || throw(ArgumentError("weight vector must be non-empty"))
   v[1] != 0 || throw(ArgumentError("first weight must be nonzero"))
   return sordering([sorder_block(t, len, v)])
end

@doc raw"""
    ordering_lp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
lexicographical ordering (:lex).
"""
ordering_lp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_lp, nvars)

@doc raw"""
    ordering_rp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
reverse lexicographical ordering (:revlex).
"""
ordering_rp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_rp, nvars)

@doc raw"""
    ordering_dp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
degree reverse lexicographical ordering (:degrevlex).
"""
ordering_dp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_dp, nvars)

@doc raw"""
    ordering_Dp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
degree lexicographical ordering (:deglex).
"""
ordering_Dp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_Dp, nvars)

@doc raw"""
    ordering_wp(w::Vector{Int})

Represents a block of variables with the
weighted reverse lexicographical ordering.
The weight vector `w` is expected to consist of positive integers only.
"""
ordering_wp(w::Vector{Int}) = _global_weighted_ordering(Singular.ringorder_wp, w)

@doc raw"""
    ordering_Wp(w::Vector{Int})

Represents a block of variables with the
weighted lexicographical ordering.
The weight vector is expected to consist of positive integers only.
"""
ordering_Wp(w::Vector{Int}) = _global_weighted_ordering(Singular.ringorder_Wp, w)

@doc raw"""
    ordering_ls(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative lexicographical ordering (:neglex).
"""
ordering_ls(nvars::Int = 1) = _basic_ordering(Singular.ringorder_ls, nvars)

@doc raw"""
    ordering_rs(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative reverse lexicographical ordering (:negrevlex).
"""
ordering_rs(nvars::Int = 1) = _basic_ordering(Singular.ringorder_rs, nvars)

@doc raw"""
    ordering_ds(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative degree reverse lexicographical ordering (:negdegrevlex).
"""
ordering_ds(nvars::Int = 1) = _basic_ordering(Singular.ringorder_ds, nvars)

@doc raw"""
    ordering_Ds(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative degree reverse lexicographical ordering (:negdeglex).
"""
ordering_Ds(nvars::Int = 1) = _basic_ordering(Singular.ringorder_Ds, nvars)

@doc raw"""
    ordering_ws(w::Vector{Int})

Represents a block of variables with the
general weighted reverse lexicographical ordering.
The weight vector `w` is expected to have a nonzero first entry.
"""
ordering_ws(w::Vector{Int}) = _local_weighted_ordering(Singular.ringorder_ws, w)

@doc raw"""
    ordering_Ws(w::Vector{Int})

Represents a block of variables with the
general weighted lexicographical ordering.
The weight vector `w` is expected to have a nonzero first entry.
"""
ordering_Ws(w::Vector{Int}) = _local_weighted_ordering(Singular.ringorder_Ws, w)

@doc raw"""
    ordering_a(w::Vector{Int})

Represents an extra weight vector that may precede any monomial ordering.
An extra weight vector does not define a monomial ordering by itself: it can
only be used in combination with other orderings to insert an extra line of
weights into the ordering matrix.
"""
ordering_a(w::Vector{Int}) = sordering([sorder_block(ringorder_a, 0, w)])

@doc raw"""
    ordering_M(m::Matrix{Int}; checked::Bool = true)

Represents a block of variables with a general matrix ordering.
The matrix `m` is expected to be invertible, and this is checked by default.
"""
function ordering_M(m::Matrix{Int}; check::Bool=true)
   (nr, nc) = size(m)
   nr > 0 && nr == nc || throw(ArgumentError("weight matrix must be square"))
   !check || !iszero(Nemo.det(Nemo.matrix(Nemo.ZZ, m))) || throw(ArgumentError("weight matrix must nonsingular"))
   return sordering([sorder_block(ringorder_M, nr, vec(transpose(m)))])
end

function ordering_M(m::ZZMatrix; check::Bool=true)
   !check || !iszero(Nemo.det(m)) || throw(ArgumentError("weight matrix must nonsingular"))
   return ordering_M(Int.(m); check=false)
end

# C, c, and S can take a dummy int in singular, but they do nothing with it?

@doc raw"""
    ordering_C()

Represents an ascending ordering on vector components `gen(1) < gen(2) < ...`.
All monomial block orderings preceding the component ordering have higher
precedence, and all succeeding monomial block orderings have lower precedence.
It is not necessary to specify this ordering explicitly since it appended
automatically to an ordering lacking a component specification.
"""
ordering_C(dummy::Int = 0) = _basic_ordering(Singular.ringorder_C, 0)

@doc raw"""
    ordering_c()

Represents a descending ordering on vector components `gen(1) > gen(2) > ...`.
All monomial block orderings preceding the component ordering have higher
precedence, and all succeeding monomial block orderings have lower precedence.
"""
ordering_c(dummy::Int = 0) = _basic_ordering(Singular.ringorder_c, 0)

ordering_S(dummy::Int = 0) = _basic_ordering(Singular.ringorder_S, 0)
ordering_s(syz_comp::Int = 0) = sordering([sorder_block(Singular.ringorder_s, syz_comp, Int[])])

@doc raw"""
    *(a::sordering, b::sordering)

Return the concatenation two orderings. Some simplification may take place,
i.e. ordering_lp(2)*ordering_lp(3) may return ordering_lp(5)
"""
function Base.:*(a::sordering, b::sordering)
   return sordering(vcat(a.data, b.data))
end

function _ispure_block(a::sordering)
   if length(a.data) == 1
      return true
   elseif length(a.data) == 2
      return a.data[2].order == ringorder_C
   else
      return false
   end
end

@doc raw"""
    ordering_size(a::sordering)

Return the size of the block of the ordering `a`, which must be a pure block.
"""
function ordering_size(a::sordering)
   _ispure_block(a) || error("ordering must be a pure block")
   return a.data[1].size
end

@doc raw"""
    ordering_weights(a::sordering)

Return the weights of the ordering `a`, which must be a pure block.
Note that for a block with an ordering specified by a matrix,
`ordering_as_symbol(a)` will return `:matrix` and the return of
`ordering_weights(a)` can be reshaped into a square matrix of dimension
`ordering_size(a)`.
"""
function ordering_weights(a::sordering)
   _ispure_block(a) || error("ordering must be a pure block")
   return a.data[1].weights
end

is_ordering_symbolic(a::sordering) = is_ordering_symbolic_with_symbol(a)[1]

@doc raw"""
    ordering_as_symbol(a::sordering)

If the ordering `a` is a pure block, return a symbol representing its type.
The symbol `:unknown` is returned if `a` is not a pure block.
"""
ordering_as_symbol(a::sordering) = is_ordering_symbolic_with_symbol(a)[2]

function is_ordering_symbolic_with_symbol(a::sordering)
   _ispure_block(a) || return (false, :unknown)
   o = a.data[1].order
   if o == ringorder_lp
      return (true, :lex)
   elseif o == ringorder_rp
      return (true, :revlex)
   elseif o == ringorder_ls
      return (true, :neglex)
   elseif o == ringorder_rs
      return (true, :negrevlex)
   elseif o == ringorder_dp
      return (true, :degrevlex)
   elseif o == ringorder_Dp
      return (true, :deglex)
   elseif o == ringorder_ds
      return (true, :negdegrevlex)
   elseif o == ringorder_Ds
      return (true, :negdeglex)
   elseif o == ringorder_wp
      return (true, :wdegrevlex)
   elseif o == ringorder_Wp
      return (true, :wdeglex)
   elseif o == ringorder_ws
      return (false, :negwdegrevlex)
   elseif o == ringorder_Ws
      return (false, :negwdeglex)
   elseif o == ringorder_a
      return (false, :extraweight)
   elseif o == ringorder_M
      return (false, :matrix)
   elseif o == ringorder_c
      return (false, :comp1max)
   elseif o == ringorder_C
      return (false, :comp1min)
   else
      return (false, :unknown)
   end
end

function Base.eltype(a::sordering)
   return sordering
end

function Base.length(a::sordering)
   return length(a.data)
end

function Base.iterate(a::sordering, state = 1)
   if state > length(a.data)
      return nothing
   else
      return sordering([a.data[state]]), state + 1
   end
end

function serialize_ordering(nvars::Int, ord::sordering)
   b = Cint[length(ord.data)]
   lastvar = 0
   cC_count = 0
   for l in 1:length(ord.data)
      i = ord.data[l]
      if i.order == ringorder_c || i.order == ringorder_C
         cC_count += 1
         if cC_count == 1
            push!(b, libSingular.ringorder_to_int(i.order))
            push!(b, 0)
            push!(b, 0)
            push!(b, 0)
         else
            error("more than one ordering c/C specified")
         end
      elseif i.order == ringorder_s || i.order == ringorder_S
         push!(b, libSingular.ringorder_to_int(i.order))
         push!(b, i.size) # blk0 and
         push!(b, i.size) # blk1 set to syz_comp for ringorder_s
         push!(b, 0)
      elseif i.order == ringorder_a
         push!(b, libSingular.ringorder_to_int(i.order))
         push!(b, lastvar + 1)
         nweights = min(length(i.weights), nvars - lastvar)
         push!(b, lastvar + nweights)
         push!(b, length(i.weights))
         for j in 1:nweights
            push!(b, i.weights[j])
         end
      else
         blksize = i.size
         if _is_weighted_ordering(i.order)
            @assert blksize > 0 && length(i.weights) == blksize
         elseif i.order == ringorder_M
            @assert blksize > 0 && length(i.weights) == blksize*blksize
         elseif _is_basic_ordering(i.order)
            @assert length(i.weights) == 0
            @assert blksize >= 0
            # consume all remaining variables when succeeded by only C,c,S,s,IS
            at_end = true
            for ll in l+1:length(ord.data)
               o = ord.data[ll].order
               if o != ringorder_C && o != ringorder_c && o != ringorder_S &&
                                      o != ringorder_s && o != ringorder_IS
                  at_end = false
                  break
               end
            end
            if at_end
               blksize = max(blksize, nvars - lastvar)
            end
         else
            error("unknown ordering $(i.order)")
         end
         push!(b, libSingular.ringorder_to_int(i.order))
         push!(b, lastvar + 1)
         lastvar += blksize
         push!(b, lastvar)
         push!(b, length(i.weights))
         for j in i.weights
            push!(b, j)
         end
      end
   end

   if nvars != lastvar
      error("mismatch of number of variables (", nvars, ") and ordering (", lastvar, ")")
   end

   # add order C if none exists
   if cC_count == 0
      b[1] += 1
      push!(b, libSingular.ringorder_to_int(ringorder_C))
      push!(b, 0)
      push!(b, 0)
      push!(b, 0)
   end

   return b
end

function deserialize_ordering(b::Vector{Cint})
   off = 0
   nblocks = b[off+=1]
   data = sorder_block[]
   for i in 1:nblocks
      o = libSingular.ringorder_from_int(b[off+=1])
      blk0 = b[off+=1]
      blk1 = b[off+=1]
      nweights = b[off+=1]
      weights = Int.(b[(off+1):(off+nweights)])
      off += nweights
      blksize = blk1 - blk0 + 1
      if _is_basic_ordering(o)
         @assert nweights == 0
         push!(data, sorder_block(o, blksize, Int[]))
      elseif _is_weighted_ordering(o) || o == ringorder_M || o == ringorder_a
         @assert nweights > 0
         @assert nweights == (o == ringorder_M ? blksize*blksize : blksize)
         push!(data, sorder_block(o, blksize, weights))
      elseif o == ringorder_C || o == ringorder_c || o == ringorder_S
         @assert nweights == 0
         push!(data, sorder_block(o, 0, Int[]))
      elseif o == ringorder_s
         push!(data, sorder_block(o, blk0, Int[]))
      else
         error("unknown ordering $o")
      end
    end
    return sordering(data)
end


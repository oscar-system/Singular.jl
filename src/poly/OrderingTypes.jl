###############################################################################
#
# sordering - for fancy orderings
#
###############################################################################

# instead of many types for each block type,
# there is one type with possibly meaningless entries
struct sorder_block
   order::Singular.libSingular.rRingOrder_t
   size::Int            # per-order meaning.
   weights::Vector{Int} # per-order meaning. empty for dp, Dp, ...
end

struct sordering
   data::Vector{sorder_block}
end


module Singular

import AbstractAlgebra
using Markdown
using Nemo

import Base: abs, checkbounds, convert, deepcopy, deepcopy_internal,
             denominator, div, divrem, exponent,
             gcd, gcdx, getindex, hash, inv, isequal, isless, lcm,
             length, mod, numerator, one, reduce, rem, setindex!, show,
             zero, +, -, *, ==, ^, &, |, <<, >>, ~, <=, >=, <, >, //,
             /, !=

import LinearAlgebra: normalize!, rank

import Statistics: std

import Nemo: add!, addeq!, base_ring, canonical_unit,
             change_base_ring, characteristic, check_parent, codomain,
             coeff, coeffs, compose, contains, content, crt,
             deflate, deflation, degree, degrees, derivative, divexact,
             divides, domain, elem_type, evaluate, exponent_vectors, finish,
             gcdinv, gen, gens, get_field, intersect, isconstant,
             isgen, ismonomial, inflate, isnegative, isone,
             isterm, isunit, iszero, lc, lt, lm, monomials,
             MPolyBuildCtx, mul!, needs_parentheses,
             nvars, ordering, parent_type,
             parent, primpart, promote_rule, push_term!,
             reconstruct, remove, set_field!, show_minus_one, sort_terms!,
             symbols, terms, total_degree, valuation,
             var_index, vars, zero!, ResidueRing

export base_ring, elem_type, parent_type, parent

export ResidueRing, PolynomialRing, Ideal, MaximalIdeal, FreeModule

export ZZ, QQ, FiniteField, FunctionField, CoefficientRing, Fp

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

const pkgdir = realpath(joinpath(dirname(@__FILE__), ".."))
const libsingular = joinpath(pkgdir, "local", "lib", "libSingular")

prefix = realpath(joinpath(@__DIR__, "..", "local"))

mapping_types = nothing
mapping_types_reversed = nothing

function __init__()

   # Initialise Singular

   binSingular = joinpath(prefix, "bin", "Singular")
   ENV["SINGULAR_EXECUTABLE"] = binSingular
   libSingular.siInit(binSingular)
   # set up Singular parents (we cannot do this before Singular is initialised)

   global ZZ = Integers()
   global QQ = Rationals()

   # done in __init__ since headers must be included first

   global n_Z_2_n_Q = libSingular.n_SetMap(ZZ.ptr, QQ.ptr)
   global n_Q_2_n_Z = libSingular.n_SetMap(QQ.ptr, ZZ.ptr)
   ZZ.refcount += 1
   QQ.refcount += 1

   global ringorder_no = libSingular.ringorder_no
   global ringorder_lp = libSingular.ringorder_lp
   global ringorder_rp = libSingular.ringorder_rp
   global ringorder_dp = libSingular.ringorder_dp
   global ringorder_Dp = libSingular.ringorder_Dp
   global ringorder_ls = libSingular.ringorder_ls
   global ringorder_rs = libSingular.ringorder_rs
   global ringorder_ds = libSingular.ringorder_ds
   global ringorder_Ds = libSingular.ringorder_Ds
   global ringorder_c  = libSingular.ringorder_c
   global ringorder_C  = libSingular.ringorder_C

   global sym2ringorder = Dict{Symbol, libSingular.rRingOrder_t}(
     :lex => ringorder_lp,
     :revlex => ringorder_rp,
     :neglex => ringorder_ls,
     :negrevlex => ringorder_rs,
     :degrevlex => ringorder_dp,
     :deglex => ringorder_Dp,
     :negdegrevlex => ringorder_ds,
     :negdeglex => ringorder_Ds,
     :comp1max => ringorder_c,
     :comp1min => ringorder_C
   )
   global mapping_types, mapping_types_reversed, casting_functions
   mapping_types = Dict( i[1] => i[2] for i in libSingular.get_type_mapper() )
   mapping_types_reversed = Dict( j => i for (i, j) in mapping_types )
   casting_functions = create_casting_functions()
end

###############################################################################
#
#   Load Singular Rings/Fields/etc
#
###############################################################################

include("AbstractTypes.jl")

include("LibSingular.jl")

include("Number.jl")

include("Poly.jl")

include("Module.jl")

include("Ideal.jl")

include("Matrix.jl")

include("Vector.jl")

include("Resolution.jl")

include("caller.jl")

include("Meta.jl")

include("Map.jl")

end # module

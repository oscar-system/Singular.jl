module Singular

import AbstractAlgebra
using Nemo

using CxxWrap

import Base: abs, checkbounds, deepcopy, deepcopy_internal,
             denominator, div, divrem, exponent,
             gcd, gcdx, getindex, indices, inv, isequal, isless, lcm, length,
             mod, normalize!, numerator, one, rank, reduce, rem, setindex!, show,
             std, zero, +, -, *, ==, ^, &, |, $, <<, >>, ~, <=, >=, <, >, //,
             /, !=

import Nemo: add!, addeq!, base_ring, canonical_unit, check_parent, coeff,
             contains, content, crt, divexact,
             elem_type, gcdinv, gens, isnegative, isone,
             isgen, iszero, isunit, lead, mul!, needs_parentheses, parent_type,
             parent, primpart, promote_rule, reconstruct, show_minus_one, zero!,
             ResidueRing, characteristic, degree, intersect, isconstant

export base_ring, elem_type, parent_type, parent

export ResidueRing, PolynomialRing, Ideal, MaximalIdeal, FreeModule

export ZZ, QQ, FiniteField, CoefficientRing, Fp

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

const pkgdir = realpath(joinpath(dirname(@__FILE__), ".."))
const libsingular = joinpath(pkgdir, "local", "lib", "libSingular")

prefix = joinpath(Pkg.dir("Singular"), "local");
nemoinc = joinpath(Pkg.dir("Nemo"), "local");

function __init__()

   # Initialise Singular
   
   binSingular = joinpath(prefix, "bin", "Singular")
   ENV["SINGULAR_EXECUTABLE"] = binSingular
   libSingular.siInit(binSingular)
   # set up Singular parents (we cannot do this before Singular is initialised)

   ZZ.ptr = get_n_Z()
   ZZ.refcount = 1

   QQ.ptr = get_n_Q()
   QQ.refcount = 1

   # done in __init__ since headers must be included first

#   global const n_Z_2_n_Q = libSingular.n_SetMap(ZZ.ptr, QQ.ptr)
#   global const n_Q_2_n_Z = libSingular.n_SetMap(QQ.ptr, ZZ.ptr)

  global const ringorder_no = libSingular.ringorder_no
  global const ringorder_lp = libSingular.ringorder_lp
  global const ringorder_rp = libSingular.ringorder_rp
  global const ringorder_dp = libSingular.ringorder_dp
  global const ringorder_Dp = libSingular.ringorder_Dp
  global const ringorder_ls = libSingular.ringorder_ls
  global const ringorder_rs = libSingular.ringorder_rs
  global const ringorder_ds = libSingular.ringorder_ds
  global const ringorder_Ds = libSingular.ringorder_Ds
  global const ringorder_c  = libSingular.ringorder_c
  global const ringorder_C  = libSingular.ringorder_C

 global const sym2ringorder = Dict{Symbol, libSingular.rRingOrder_t}(
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
end

###############################################################################
#
#   Load Singular Rings/Fields/etc
#
###############################################################################

include("AbstractTypes.jl")

include("LibSingular.jl")

include("Number.jl")

# include("Poly.jl")

# include("Module.jl")

# include("Ideal.jl")

# include("Matrix.jl")

# include("Vector.jl")

# include("Resolution.jl")

###############################################################################
#
#   Set Singlular ZZ and QQ
#
###############################################################################

ZZ = Integers()

QQ = Rationals()

end # module

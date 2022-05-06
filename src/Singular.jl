@doc Markdown.doc"""
  Julia package for using the Singular library for commutative and non-commutative algebra, algebraic geometry, and singularity theory.

  For documentation see https://oscar-system.github.io/Singular.jl/latest/index.html

  For more information about Singular see https://www.singular.uni-kl.de/
"""
module Singular

#
# More helpful error message for users on Windows.
windows_error() = error("""

    This package unfortunately does not run natively under Windows.
    Please install Julia using Windows subsystem for Linux and try again.
    See also https://oscar.computeralgebra.de/install/.
    """)

if Sys.iswindows()
  windows_error()
end


import AbstractAlgebra
using Markdown
using Nemo
using Pkg
using lib4ti2_jll
using Singular_jll
using libsingular_julia_jll

import Base: abs, checkbounds, convert, deepcopy, deepcopy_internal,
             denominator, div, divrem, exponent,
             gcd, gcdx, getindex, hash, inv, isequal, isless, lcm,
             length, mod, numerator, one, reduce, rem, setindex!, show,
             transpose, zero, +, -, *, ==, ^, &, |, <<, >>, ~, <=,
             >=, <, >, //, /, !=

import LinearAlgebra: normalize!, rank

import Statistics: std

using Random: Random, AbstractRNG, SamplerTrivial, SamplerSimple
import Random: rand
using RandomExtensions: RandomExtensions, make, Make2

import AbstractAlgebra: AbstractAlgebra, diagonal_matrix, factor,
                        identity_matrix, kernel, lead, ncols, ngens, nrows, order,
                        preimage, zero_matrix, expressify

import Nemo: add!, addeq!, base_ring, canonical_unit,
             change_base_ring, characteristic, check_parent, codomain,
             coeff, coefficients, compose, constant_coefficient,
             contains, content, crt,
             deflate, deflation, degree, degrees, derivative, divexact,
             divides, domain, elem_type, evaluate, exponent_vectors,
             factor_squarefree,
             finish, gcdinv, gen, gens, intersect, is_constant,
             is_gen, is_monomial, inflate, is_negative, isone,
             is_term, is_unit, iszero, lift, leading_coefficient,
             leading_term, leading_monomial, monomials,
             MPolyBuildCtx, mul!, nvars, ordering, parent_type,
             parent, primpart, promote_rule, push_term!,
             reconstruct, remove, sort_terms!,
             symbols, tail, terms, total_degree, trailing_coefficient, valuation,
             var_index, vars, zero!, ResidueRing

export base_ring, elem_type, parent_type, parent

export ResidueRing, PolynomialRing, ExteriorAlgebra, WeylAlgebra, Ideal,
       MaximalIdeal, FreeModule, @PolynomialRing, @WeylAlgebra,
       @ExteriorAlgebra

export ZZ, QQ, FiniteField, FunctionField, CoefficientRing, Fp

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

const libsingular = Singular_jll.libsingular
const binSingular = Singular_jll.Singular_path

const libsingula_julia = libsingular_julia_jll.libsingular_julia

const libflint = Nemo.libflint
const libantic = Nemo.libantic

mapping_types_reversed = nothing

function __init__()
   if Sys.iswindows()
      windows_error()
   end

   # Initialise Singular
   ENV["SINGULAR_EXECUTABLE"] = binSingular
   libSingular.siInit(binSingular)

   # don't tell the user about each library that singular loads
   libSingular.set_option("V_LOAD_LIB", false)
   # don't tell the user when a library redefines a singular variable
   libSingular.set_option("V_REDEFINE", false)
   # silence printlevel-based printing
   libSingular.set_printlevel(-100)

   # At this point we have Kstd1_mu = 2^31-1 and OPT_MULTBOUND unset, which is
   # a slightly inconsistent state. The singular intepreter variable "multBound"
   # simply reads Kstd1_mu and assigning to "multBound" assigns to Kstd1_mu
   # and sets OPT_MULTBOUND if it is nonzero, and unsets OPT_MULTBOUND if it is
   # zero. Since we want to be able to easily restore the multBound to a
   # consistent state, we set Kstd1_mu to zero here, which also leaves
   # OPT_MULTBOUND unset.
   Singular.libSingular.set_multBound(0)
   # ditto
   Singular.libSingular.set_degBound(0)

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
   global ringorder_wp = libSingular.ringorder_wp
   global ringorder_Wp = libSingular.ringorder_Wp
   global ringorder_ls = libSingular.ringorder_ls
   global ringorder_rs = libSingular.ringorder_rs
   global ringorder_ds = libSingular.ringorder_ds
   global ringorder_Ds = libSingular.ringorder_Ds
   global ringorder_ws = libSingular.ringorder_ws
   global ringorder_Ws = libSingular.ringorder_Ws
   global ringorder_a  = libSingular.ringorder_a
   global ringorder_M  = libSingular.ringorder_M
   global ringorder_c  = libSingular.ringorder_c
   global ringorder_C  = libSingular.ringorder_C
   global ringorder_s  = libSingular.ringorder_s
   global ringorder_S  = libSingular.ringorder_S
   global ringorder_IS = libSingular.ringorder_IS

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
   global mapping_types_reversed, casting_functions
   mapping_types_reversed = Dict( i[2] => i[1] for i in libSingular.get_type_mapper() )
   casting_functions = create_casting_functions()

   # Respect the -q and --banner flag
   allowbanner = Base.JLOptions().banner != 0

   show_banner = allowbanner && isinteractive() &&
                !any(x->x.name in ["Oscar"], keys(Base.package_locks)) &&
                get(ENV, "SINGULAR_PRINT_BANNER", "true") != "false"

   singular_version_nr=Singular.libSingular.version()
   ver = digits(singular_version_nr, base = 10)
   svn = "$(ver[4]).$(ver[3]).$(ver[2])"
   if ver[1] > 0
     svn *= "p$(ver[1])"
   end
   if show_banner
     println("""Singular.jl, based on
                     SINGULAR                               /
 A Computer Algebra System for Polynomial Computations     /  Singular.jl: $VERSION_NUMBER
                                                         0<   Singular   : $svn
 by: W. Decker, G.-M. Greuel, G. Pfister, H. Schoenemann   \\
FB Mathematik der Universitaet, D-67653 Kaiserslautern      \\
     """)
   end
end

pkgproject(m::Core.Module) = Pkg.Operations.read_project(Pkg.Types.projectfile_path(pkgdir(m)))
pkgversion(m::Core.Module) = pkgproject(m).version
const VERSION_NUMBER = pkgversion(@__MODULE__)


###############################################################################
#
#   Load Singular Rings/Fields/etc
#     There is a slight circular dependency: in order to create a G-Algebra
#     you need a PolyRing R and matrices over R, and in order to create a
#     quotient ring you need the polynomial ring and an ideal.
#     The matrix code requires some module code, which requires some ideal
#     code, which requies some poly code. Therefore, we include all of the
#     poly/matrix/module/ideal types before including the methods.
#
###############################################################################

include("setup.jl")

include("AbstractTypes.jl")

include("LibSingular.jl")

include("Number.jl")

include("poly/OrderingTypes.jl")

include("poly/PolyTypes.jl")

include("matrix/MatrixTypes.jl")

include("module/ModuleTypes.jl")

include("poly/PluralTypes.jl")

include("poly/LPTypes.jl")

# all "poly" types

const PolyRingUnion{T} = Union{PolyRing{T}, PluralRing{T}, LPRing{T}} where T <: Nemo.RingElem

const SPolyUnion{T} = Union{spoly{T}, spluralg{T}, slpalg{T}} where T <: Nemo.RingElem

include("poly/orderings.jl")

include("poly/poly.jl")

include("poly/plural.jl")

include("poly/weyl.jl")

include("poly/lp.jl")

include("ideal/IdealTypes.jl")

include("matrix/matrix.jl")

include("matrix/bigintmat.jl")

include("ideal/quotient.jl")

include("module/module.jl")

include("ideal/ideal.jl")

include("Vector.jl")

include("Resolution.jl")

include("caller.jl")

include("Meta.jl")

include("Map.jl")

include("MessyHacks.jl")

end # module
